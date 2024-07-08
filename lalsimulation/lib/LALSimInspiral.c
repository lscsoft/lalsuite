/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, D. Keppel, Evan Ochsner, C. Pankow, 2014 A. Klein
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

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralEOS.h>
#include <lal/LALSimIMR.h>
#include <lal/LALSimSphHarmMode.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Units.h>
#include <lal/LALSimBlackHoleRingdown.h>
#include <lal/LALSimInspiralPrecess.h>
#include <lal/LALSimInspiralWaveformParams.h>

#include "LALSimInspiralPNCoefficients.c"
#include "check_series_macros.h"
#include "check_waveform_macros.h"
#include "LALSimUniversalRelations.h"
#include "rotation_macros.h"
#include "fix_reference_frequency_macro.h"

#include "LALSimInspiralGenerator_private.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * (Twice) the highest known PN order of amplitude correction for
 * non-precessing binaries.
 */
#define MAX_NONPRECESSING_AMP_PN_ORDER 6

/**
 * (Twice) the highest known PN order of amplitude correction for
 * precessing binaries.
 */
#define MAX_PRECESSING_AMP_PN_ORDER 3




#define INITIALIZE_NAME(a) [a] = #a
/* TODO: UPDATE ME WHENEVER A NEW APPROXIMANT IS ADDED */
static const char *lalSimulationApproximantNames[] = {
    INITIALIZE_NAME(TaylorT1),
    INITIALIZE_NAME(TaylorT2),
    INITIALIZE_NAME(TaylorT3),
    INITIALIZE_NAME(TaylorF1),
    INITIALIZE_NAME(TaylorF2),
    INITIALIZE_NAME(TaylorF2Ecc),
    INITIALIZE_NAME(TaylorF2NLTides),
    INITIALIZE_NAME(TaylorR2F4),
    INITIALIZE_NAME(TaylorF2RedSpin),
    INITIALIZE_NAME(TaylorF2RedSpinTidal),
    INITIALIZE_NAME(PadeT1),
    INITIALIZE_NAME(PadeF1),
    INITIALIZE_NAME(EOB),
    INITIALIZE_NAME(BCV),
    INITIALIZE_NAME(BCVSpin),
    INITIALIZE_NAME(SpinTaylorT1),
    INITIALIZE_NAME(SpinTaylorT5),
    INITIALIZE_NAME(SpinTaylorT3),
    INITIALIZE_NAME(SpinTaylorT4),
    INITIALIZE_NAME(SpinTaylorF2),
    INITIALIZE_NAME(SpinTaylorFrameless),
    INITIALIZE_NAME(SpinTaylor),
    INITIALIZE_NAME(PhenSpinTaylor),
    INITIALIZE_NAME(PhenSpinTaylorRD),
    INITIALIZE_NAME(SpinQuadTaylor),
    INITIALIZE_NAME(FindChirpSP),
    INITIALIZE_NAME(FindChirpPTF),
    INITIALIZE_NAME(GeneratePPN),
    INITIALIZE_NAME(BCVC),
    INITIALIZE_NAME(FrameFile),
    INITIALIZE_NAME(AmpCorPPN),
    INITIALIZE_NAME(NumRel),
    INITIALIZE_NAME(NumRelNinja2),
    INITIALIZE_NAME(EccentricFD),
    INITIALIZE_NAME(Eccentricity),
    INITIALIZE_NAME(EOBNR),
    INITIALIZE_NAME(EOBNRv2),
    INITIALIZE_NAME(EOBNRv2HM),
    INITIALIZE_NAME(EOBNRv2_ROM),
    INITIALIZE_NAME(EOBNRv2HM_ROM),
    INITIALIZE_NAME(TEOBResum_ROM),
    INITIALIZE_NAME(SEOBNRv1),
    INITIALIZE_NAME(SEOBNRv2),
    INITIALIZE_NAME(SEOBNRv2_opt),
    INITIALIZE_NAME(SEOBNRv3),
    INITIALIZE_NAME(SEOBNRv3_pert),
    INITIALIZE_NAME(SEOBNRv3_opt),
    INITIALIZE_NAME(SEOBNRv3_opt_rk4),
    INITIALIZE_NAME(SEOBNRv4),
    INITIALIZE_NAME(SEOBNRv4_opt),
    INITIALIZE_NAME(SEOBNRv4P),
    INITIALIZE_NAME(SEOBNRv4PHM),
    INITIALIZE_NAME(SEOBNRv2T),
    INITIALIZE_NAME(SEOBNRv4T),
    INITIALIZE_NAME(SEOBNRv4HM),
    INITIALIZE_NAME(SEOBNRv1_ROM_EffectiveSpin),
    INITIALIZE_NAME(SEOBNRv1_ROM_DoubleSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_EffectiveSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_DoubleSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_DoubleSpin_HI),
    INITIALIZE_NAME(Lackey_Tidal_2013_SEOBNRv2_ROM),
    INITIALIZE_NAME(SEOBNRv4_ROM),
    INITIALIZE_NAME(SEOBNRv4HM_ROM),
    INITIALIZE_NAME(SEOBNRv4_ROM_NRTidal),
    INITIALIZE_NAME(SEOBNRv4_ROM_NRTidalv2),
    INITIALIZE_NAME(SEOBNRv4_ROM_NRTidalv2_NSBH),
    INITIALIZE_NAME(SEOBNRv4T_surrogate),
    INITIALIZE_NAME(SEOBNRv5_ROM),
    INITIALIZE_NAME(SEOBNRv5HM_ROM),
    INITIALIZE_NAME(SEOBNRv5_ROM_NRTidalv3),
    INITIALIZE_NAME(HGimri),
    INITIALIZE_NAME(IMRPhenomA),
    INITIALIZE_NAME(IMRPhenomB),
    INITIALIZE_NAME(IMRPhenomFA),
    INITIALIZE_NAME(IMRPhenomFB),
    INITIALIZE_NAME(IMRPhenomC),
    INITIALIZE_NAME(IMRPhenomD),
    INITIALIZE_NAME(IMRPhenomD_NRTidal),
    INITIALIZE_NAME(IMRPhenomD_NRTidalv2),
    INITIALIZE_NAME(IMRPhenomNSBH),
    INITIALIZE_NAME(IMRPhenomHM),
    INITIALIZE_NAME(IMRPhenomP),
    INITIALIZE_NAME(IMRPhenomPv2),
    INITIALIZE_NAME(IMRPhenomPv2_NRTidal),
    INITIALIZE_NAME(IMRPhenomPv2_NRTidalv2),
    INITIALIZE_NAME(IMRPhenomPv3),
    INITIALIZE_NAME(IMRPhenomPv3HM),
    INITIALIZE_NAME(IMRPhenomFC),
    INITIALIZE_NAME(TaylorEt),
    INITIALIZE_NAME(TaylorT4),
    INITIALIZE_NAME(EccentricTD),
    INITIALIZE_NAME(TaylorN),
    INITIALIZE_NAME(SpinTaylorT4Fourier),
    INITIALIZE_NAME(SpinTaylorT5Fourier),
    INITIALIZE_NAME(SpinDominatedWf),
    INITIALIZE_NAME(NRSur4d2s),
    INITIALIZE_NAME(NRSur7dq2),
    INITIALIZE_NAME(NRSur7dq4),
    INITIALIZE_NAME(NR_hdf5),
    INITIALIZE_NAME(NRHybSur3dq8),
    INITIALIZE_NAME(IMRPhenomXAS),
    INITIALIZE_NAME(IMRPhenomXHM),
    INITIALIZE_NAME(IMRPhenomXP),
    INITIALIZE_NAME(IMRPhenomXPHM),
    INITIALIZE_NAME(TEOBResumS),
    INITIALIZE_NAME(IMRPhenomT),
    INITIALIZE_NAME(IMRPhenomTHM),
    INITIALIZE_NAME(IMRPhenomTP),
    INITIALIZE_NAME(IMRPhenomTPHM),
    INITIALIZE_NAME(IMRPhenomXO4a),
    INITIALIZE_NAME(IMRPhenomXAS_NRTidalv2),
    INITIALIZE_NAME(IMRPhenomXP_NRTidalv2),
    INITIALIZE_NAME(IMRPhenomXAS_NRTidalv3),
    INITIALIZE_NAME(IMRPhenomXP_NRTidalv3),
    INITIALIZE_NAME(SEOBNRv4HM_PA),
    INITIALIZE_NAME(pSEOBNRv4HM_PA),
    INITIALIZE_NAME(ExternalPython)
};
#undef INITIALIZE_NAME

/* TODO: UPDATE ME WHENEVER A NEW PN ORDER IS ADDED */
static const char *lalSimulationPNOrderNames[] = {
    [LAL_PNORDER_NEWTONIAN]         = "newtonian",
    [LAL_PNORDER_HALF]              = "oneHalfPN",
    [LAL_PNORDER_ONE]               = "onePN",
    [LAL_PNORDER_ONE_POINT_FIVE]    = "onePointFivePN",
    [LAL_PNORDER_TWO]               = "twoPN",
    [LAL_PNORDER_TWO_POINT_FIVE]    = "twoPointFivePN",
    [LAL_PNORDER_THREE]             = "threePN",
    [LAL_PNORDER_THREE_POINT_FIVE]  = "threePointFivePN",
    [LAL_PNORDER_PSEUDO_FOUR]       = "pseudoFourPN",
};

/* TODO: UPDATE ME WHENEVER A NEW TAPER IS ADDED */
static const char *lalSimulationTaperNames[] = {
    [LAL_SIM_INSPIRAL_TAPER_NONE]       = "TAPER_NONE",
    [LAL_SIM_INSPIRAL_TAPER_START]      = "TAPER_START",
    [LAL_SIM_INSPIRAL_TAPER_END]        = "TAPER_END",
    [LAL_SIM_INSPIRAL_TAPER_STARTEND]   = "TAPER_STARTEND",
};

/* TODO: UPDATE ME WHENEVER A NEW FRAME AXIS IS ADDED */
static const char *lalSimulationFrameAxisNames[] = {
    [LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J]   = "TotalJ",
    [LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L] = "OrbitalL",
    [LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW]      = "View",
};

/* TODO: UPDATE ME WHENEVER A NEW MODES CHOICE IS ADDED */
static const char *lalSimulationModesChoiceNames[] = {
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3AND4AND5L] = "L2345",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3AND4L] = "L234",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3AND5L] = "L235",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND4AND5L] = "L245",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_3AND4AND5L] = "L345",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3L] = "L23",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND4L] = "L24",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_3AND4L] = "L34",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_2AND5L] = "L25",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_3AND5L] = "L35",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_4AND5L] = "L45",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_RESTRICTED] = "L2",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_3L] = "L3",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_4L] = "L4",
    [LAL_SIM_INSPIRAL_MODES_CHOICE_5L] = "L5",
    /* NOTE: cannot do the "ALL" case since its value is -1 */
    // [LAL_SIM_INSPIRAL_MODES_CHOICE_ALL] = "ALL",
};

/* locates and deletes a substring in a list of substrings from a string, ignoring case;
 * if multiple substrings in the string match, delete the longest one; here, deletion
 * means replacing the substring with BEL characters */
static int delete_substring_in_list_from_string(char *string, const char *list[], size_t size)
{
    int longest_position = -1;
    int longest_offset = -1;
    int longest_length = -1;
    size_t i;

    if (string == NULL || strlen(string) == 0) // no string to search
        return -1;

    for (i = 0; i < size; ++i) {
        char *match;
        if (list[i] == NULL) // no such element in list
            continue;
        if ((match = XLALStringCaseSubstring(string, list[i]))) {
            int length = strlen(list[i]);
            if (length > longest_length) {
                longest_position = i;
                longest_offset = match - string;
                longest_length = length;
            }
        }
    }

    if (longest_position < 0) // failed to find a word
        return -1;

    /* delete word from string by replacing with BEL */
    for (i = 0; i < (size_t)longest_length; ++i)
        string[longest_offset + i] = '\b';

    return longest_position;
}

const LALSimInspiralGenerator *lalSimInspiralGeneratorTemplates[NumApproximants] = {
    [EOBNRv2HM] = &lalEOBNRv2HMGeneratorTemplate,
    [EOBNRv2HM_ROM] = &lalEOBNRv2HM_ROMGeneratorTemplate,
    [EOBNRv2] = &lalEOBNRv2GeneratorTemplate,
    [EOBNRv2_ROM] = &lalEOBNRv2_ROMGeneratorTemplate,
    [EccentricFD] = &lalEccentricFDGeneratorTemplate,
    [EccentricTD] = &lalEccentricTDGeneratorTemplate,
    [HGimri] = &lalHGimriGeneratorTemplate,
    [IMRPhenomA] = &lalIMRPhenomAGeneratorTemplate,
    [IMRPhenomB] = &lalIMRPhenomBGeneratorTemplate,
    [IMRPhenomC] = &lalIMRPhenomCGeneratorTemplate,
    [IMRPhenomD] = &lalIMRPhenomDGeneratorTemplate,
    [IMRPhenomD_NRTidal] = &lalIMRPhenomD_NRTidalGeneratorTemplate,
    [IMRPhenomD_NRTidalv2] = &lalIMRPhenomD_NRTidalv2GeneratorTemplate,
    [IMRPhenomHM] = &lalIMRPhenomHMGeneratorTemplate,
    [IMRPhenomNSBH] = &lalIMRPhenomNSBHGeneratorTemplate,
    [IMRPhenomP] = &lalIMRPhenomPGeneratorTemplate,
    [IMRPhenomPv2] = &lalIMRPhenomPv2GeneratorTemplate,
    [IMRPhenomPv2_NRTidal] = &lalIMRPhenomPv2_NRTidalGeneratorTemplate,
    [IMRPhenomPv2_NRTidalv2] = &lalIMRPhenomPv2_NRTidalv2GeneratorTemplate,
    [IMRPhenomPv3HM] = &lalIMRPhenomPv3HMGeneratorTemplate,
    [IMRPhenomPv3] = &lalIMRPhenomPv3GeneratorTemplate,
    [IMRPhenomTHM] = &lalIMRPhenomTHMGeneratorTemplate,
    [IMRPhenomTPHM] = &lalIMRPhenomTPHMGeneratorTemplate,
    [IMRPhenomTP] = &lalIMRPhenomTPGeneratorTemplate,
    [IMRPhenomT] = &lalIMRPhenomTGeneratorTemplate,
    [IMRPhenomXAS] = &lalIMRPhenomXASGeneratorTemplate,
    [IMRPhenomXHM] = &lalIMRPhenomXHMGeneratorTemplate,
    [IMRPhenomXPHM] = &lalIMRPhenomXPHMGeneratorTemplate,
    [IMRPhenomXP] = &lalIMRPhenomXPGeneratorTemplate,
    [IMRPhenomXAS_NRTidalv2] = &lalIMRPhenomXAS_NRTidalv2GeneratorTemplate,
    [IMRPhenomXP_NRTidalv2] = &lalIMRPhenomXP_NRTidalv2GeneratorTemplate,
    [IMRPhenomXAS_NRTidalv3] = &lalIMRPhenomXAS_NRTidalv3GeneratorTemplate,
    [IMRPhenomXP_NRTidalv3] = &lalIMRPhenomXP_NRTidalv3GeneratorTemplate,
    [IMRPhenomXO4a] = &lalIMRPhenomXO4aGeneratorTemplate,
    [Lackey_Tidal_2013_SEOBNRv2_ROM] = &lalLackey_Tidal_2013_SEOBNRv2_ROMGeneratorTemplate,
    [NRHybSur3dq8] = &lalNRHybSur3dq8GeneratorTemplate,
    [NRSur4d2s] = &lalNRSur4d2sGeneratorTemplate,
    [NRSur7dq2] = &lalNRSur7dq2GeneratorTemplate,
    [NRSur7dq4] = &lalNRSur7dq4GeneratorTemplate,
    [NR_hdf5] = &lalNR_hdf5GeneratorTemplate,
    [PhenSpinTaylorRD] = &lalPhenSpinTaylorRDGeneratorTemplate,
    [PhenSpinTaylor] = &lalPhenSpinTaylorGeneratorTemplate,
    [SEOBNRv1] = &lalSEOBNRv1GeneratorTemplate,
    [SEOBNRv1_ROM_DoubleSpin] = &lalSEOBNRv1_ROM_DoubleSpinGeneratorTemplate,
    [SEOBNRv1_ROM_EffectiveSpin] = &lalSEOBNRv1_ROM_EffectiveSpinGeneratorTemplate,
    [SEOBNRv2T] = &lalSEOBNRv2TGeneratorTemplate,
    [SEOBNRv2] = &lalSEOBNRv2GeneratorTemplate,
    [SEOBNRv2_ROM_DoubleSpin] = &lalSEOBNRv2_ROM_DoubleSpinGeneratorTemplate,
    [SEOBNRv2_ROM_DoubleSpin_HI] = &lalSEOBNRv2_ROM_DoubleSpin_HIGeneratorTemplate,
    [SEOBNRv2_ROM_EffectiveSpin] = &lalSEOBNRv2_ROM_EffectiveSpinGeneratorTemplate,
    [SEOBNRv2_opt] = &lalSEOBNRv2_optGeneratorTemplate,
    [SEOBNRv3] = &lalSEOBNRv3GeneratorTemplate,
    [SEOBNRv3_opt] = &lalSEOBNRv3_optGeneratorTemplate,
    [SEOBNRv3_opt_rk4] = &lalSEOBNRv3_opt_rk4GeneratorTemplate,
    [SEOBNRv3_pert] = &lalSEOBNRv3_pertGeneratorTemplate,
    [SEOBNRv4HM] = &lalSEOBNRv4HMGeneratorTemplate,
    [SEOBNRv4HM_ROM] = &lalSEOBNRv4HM_ROMGeneratorTemplate,
    [SEOBNRv4PHM] = &lalSEOBNRv4PHMGeneratorTemplate,
    [SEOBNRv4P] = &lalSEOBNRv4PGeneratorTemplate,
    [SEOBNRv4T] = &lalSEOBNRv4TGeneratorTemplate,
    [SEOBNRv4T_surrogate] = &lalSEOBNRv4T_surrogateGeneratorTemplate,
    [SEOBNRv4] = &lalSEOBNRv4GeneratorTemplate,
    [SEOBNRv4_ROM] = &lalSEOBNRv4_ROMGeneratorTemplate,
    [SEOBNRv4_ROM_NRTidal] = &lalSEOBNRv4_ROM_NRTidalGeneratorTemplate,
    [SEOBNRv4_ROM_NRTidalv2] = &lalSEOBNRv4_ROM_NRTidalv2GeneratorTemplate,
    [SEOBNRv4_ROM_NRTidalv2_NSBH] = &lalSEOBNRv4_ROM_NRTidalv2_NSBHGeneratorTemplate,
    [SEOBNRv4_opt] = &lalSEOBNRv4_optGeneratorTemplate,
    [SEOBNRv5_ROM] = &lalSEOBNRv5_ROMGeneratorTemplate,
    [SEOBNRv5HM_ROM] = &lalSEOBNRv5HM_ROMGeneratorTemplate,
    [SEOBNRv5_ROM_NRTidalv3] = &lalSEOBNRv5_ROM_NRTidalv3GeneratorTemplate,
    [SEOBNRv4HM_PA] = &lalSEOBNRv4HM_PAGeneratorTemplate,
    [pSEOBNRv4HM_PA] = &lalpSEOBNRv4HM_PAGeneratorTemplate,
    [SpinDominatedWf] = &lalSpinDominatedWfGeneratorTemplate,
    [SpinTaylorF2] = &lalSpinTaylorF2GeneratorTemplate,
    [SpinTaylorT1] = &lalSpinTaylorT1GeneratorTemplate,
    [SpinTaylorT4Fourier] = &lalSpinTaylorT4FourierGeneratorTemplate,
    [SpinTaylorT4] = &lalSpinTaylorT4GeneratorTemplate,
    [SpinTaylorT5Fourier] = &lalSpinTaylorT5FourierGeneratorTemplate,
    [SpinTaylorT5] = &lalSpinTaylorT5GeneratorTemplate,
    [TEOBResumS] = &lalTEOBResumSGeneratorTemplate,
    [TEOBResum_ROM] = &lalTEOBResum_ROMGeneratorTemplate,
    [TaylorEt] = &lalTaylorEtGeneratorTemplate,
    [TaylorF2Ecc] = &lalTaylorF2EccGeneratorTemplate,
    [TaylorF2NLTides] = &lalTaylorF2NLTidesGeneratorTemplate,
    [TaylorF2RedSpinTidal] = &lalTaylorF2RedSpinTidalGeneratorTemplate,
    [TaylorF2RedSpin] = &lalTaylorF2RedSpinGeneratorTemplate,
    [TaylorF2] = &lalTaylorF2GeneratorTemplate,
    [TaylorR2F4] = &lalTaylorR2F4GeneratorTemplate,
    [TaylorT1] = &lalTaylorT1GeneratorTemplate,
    [TaylorT2] = &lalTaylorT2GeneratorTemplate,
    [TaylorT3] = &lalTaylorT3GeneratorTemplate,
    [TaylorT4] = &lalTaylorT4GeneratorTemplate,
    [ExternalPython] = &lalPythonGeneratorTemplate,
};

/**
 * @addtogroup LALSimInspiral_c
 * @brief General routines for generating binary inspiral waveforms.
 *
 * @{
 */
 
 /**
  * @name New Interface Generator Routines
  * @{
  */
  
  /**
   * Destroy LALSimInspiralGenerator object.
   */
void XLALDestroySimInspiralGenerator(LALSimInspiralGenerator *generator)
{
    if (generator) {
        if (generator->initialize == NULL && generator->finalize == NULL) {
            /* assume generator is immutable -- this is a no-op */
            return;
        }
        /* invoke finalizer if present */
        if (generator->finalize)
            if (generator->finalize(generator) < 0)
                XLAL_ERROR_VOID(XLAL_EFUNC);
        XLALFree(generator);
    }
    return;
}

/**
 * Create LALSimInspiralGenerator object.
 * The LALDict can be None for all the C approximants. For external python model, the LALDict must contain the file and class object where the `generate_waveform` methods are defined.
 * Activating the option `condition` in the LALDict will generate the waveform with adequate conditioning for (inverse)Fourier transform.
 */
LALSimInspiralGenerator *XLALCreateSimInspiralGenerator(const LALSimInspiralGenerator *generator, LALDict *params)
{
    LALSimInspiralGenerator *new;
    XLAL_CHECK_NULL(generator, XLAL_EFAULT);
    if (generator->initialize == NULL && generator->finalize == NULL) {
        /* assume generator is immutable -- return the template */
        /* danger! discard the const qualifier! */
        new = (LALSimInspiralGenerator *)(intptr_t)generator;
    } else {
        /* create a new instance of the generator */
        new = XLALMalloc(sizeof(*new));
        XLAL_CHECK_NULL(new, XLAL_ENOMEM, "could not allocate memory for new generator");
        memcpy(new, generator, sizeof(*new));
        /* invoke initializer if present */
        if (new->initialize)
            if (new->initialize(new, params) < 0) {
                XLALFree(new);
                XLAL_ERROR_NULL(XLAL_EFUNC);
            }
    }
    return new;
}

/**
 * Returns LALSimInspiralGenerator object from approximant.
 * The LALDict can be None for all the C approximants. For external python model, the LALDict must contain the file and class object where the `generate_waveform` methods are defined.
 * Activating the option `condition` in the LALDict will generate the waveform with adequate conditioning for (inverse)Fourier transform.
 */
LALSimInspiralGenerator *XLALSimInspiralChooseGenerator(Approximant approx, LALDict *params)
{
    const LALSimInspiralGenerator *template = lalSimInspiralGeneratorTemplates[approx];
    XLAL_CHECK_NULL(template, XLAL_EINVAL, "no generator defined for approximant %d", approx);
    return XLALCreateSimInspiralGenerator(template, params);
}

/**
 * Return approximant name from generator object
 *
 */
const char *XLALSimInspiralGeneratorName(LALSimInspiralGenerator *generator)
{
    XLAL_CHECK_NULL(generator, XLAL_EFAULT);
    return generator->name;
}

 /** @} */

/**
 * @name New Interface Waveform Routines
 * @{
 */
 
/**
 * Returns time-domain polarizations for a specific approximant.
 * Equivalent to XLALSimInspiralChooseTDWaveform(). Equivalent to XLALSimInspiralTD() if the option `condition` is activated in the LALDict.
 * The waveform arguments are inserted into the LALDict. The generator carries the info about the approximant and potentially extra data which could be recycled by the model to speed-up calculation.  
 *
 * The parameters in the LALDict must be in SI units.
 */
int XLALSimInspiralGenerateTDWaveform(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *generator
)
{
    XLAL_CHECK(hplus && hcross && generator, XLAL_EFAULT);
    XLAL_CHECK(*hplus == NULL && *hcross == NULL, XLAL_EINVAL, "hplus and hcross must be pointers to NULL");

    if (generator->generate_td_waveform)
        return generator->generate_td_waveform(hplus, hcross, params, generator);

    XLAL_ERROR(XLAL_EINVAL, "generator does not provide a method to generate time-domain waveforms");
}

/**
 * Compute time-domain modes for a specific approximant.
 * Equivalent to XLALSimInspiralChooseTDModes(). The only difference is that the SphHarmSeries object needs to be passed as an argument to the function. The actual returned value is an integer which indicates success or error in the waveform evaluation (see https://lscsoft.docs.ligo.org/lalsuite/lal/group___x_l_a_l_error__h.html).
 * The waveform arguments are inserted into the LALDict. The generator carries the info about the approximant and potentially extra data which could be recycled by the model to speed-up calculation.  
 *
 * The parameters in the LALDict must be in SI units.
 */
int XLALSimInspiralGenerateTDModes(
    SphHarmTimeSeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *generator
)
{
    XLAL_CHECK(hlm && generator, XLAL_EFAULT);
    XLAL_CHECK(*hlm == NULL, XLAL_EINVAL, "hlm must be a pointer to NULL");

    if (generator->generate_td_modes)
        return generator->generate_td_modes(hlm, params, generator);

    XLAL_ERROR(XLAL_EINVAL, "generator does not provide a method to generate time-domain modes");
}

/**
 * Returns frequency-domain polarizations for a specific approximant.
 * Equivalent to XLALSimInspiralChooseFDWaveform(). Equivalent to XLALSimInspiralFD() if the option `condition` is activated in the LALDict.
 * The waveform arguments are inserted into the LALDict. The generator carries the info about the approximant and potentially extra data which could be recycled by the model to speed-up calculation.  
 *
 * The parameters in the LALDict must be in SI units.
 */
int XLALSimInspiralGenerateFDWaveform(
    COMPLEX16FrequencySeries **hplus,
    COMPLEX16FrequencySeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *generator
)
{
    XLAL_CHECK(hplus && hcross && generator, XLAL_EFAULT);
    XLAL_CHECK(*hplus == NULL && *hcross == NULL, XLAL_EINVAL, "hplus and hcross must be pointers to NULL");
    if (generator->generate_fd_waveform)
        return generator->generate_fd_waveform(hplus, hcross, params, generator);

    XLAL_ERROR(XLAL_EINVAL, "generator does not provide a method to generate frequency-domain waveforms");
}

/**
 * Compute frequency-domain modes for a specific approximant. 
 * Equivalent to XLALSimInspiralChooseFDModes. The only difference is that the SphHarmSeries object needs to be passed as an argument to the function. The actual returned value is an integer which indicates success or error in the waveform evaluation (see https://lscsoft.docs.ligo.org/lalsuite/lal/group___x_l_a_l_error__h.html).
 * The waveform arguments are inserted into the LALDict. The generator carries the info about the approximant and potentially extra data which could be recycled by the model to speed-up calculation.  
 *
 * The parameters in the LALDict must be in SI units.
 */
int XLALSimInspiralGenerateFDModes(
    SphHarmFrequencySeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *generator
)
{
    XLAL_CHECK(hlm && generator, XLAL_EFAULT);
    XLAL_CHECK(*hlm == NULL, XLAL_EINVAL, "hlm must be a pointer to NULL");

    if (generator->generate_fd_modes)
        return generator->generate_fd_modes(hlm, params, generator);

    XLAL_ERROR(XLAL_EINVAL, "generator does not provide a method to generate frequency-domain modes");
}

/** @} */

/**
* @name New Interface parse parameters Routines
* @{
*/

/**
 * Insert all the input arguments needed by XALSimInspiralChooseTDWaveform() into a laldictionary.
 */
void XLALSimInspiralParseDictionaryToChooseTDWaveform(
    REAL8 *m1,                             /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                             /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                            /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                            /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                            /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                            /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                            /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                            /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *distance,                       /**< [out] distance of source (m) */
    REAL8 *inclination,                    /**< [out] inclination of source (rad) */
    REAL8 *phiRef,                         /**< [out] reference orbital phase (rad) */
    REAL8 *longAscNodes,                   /**< [out] longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 *eccentricity,                   /**< [out] eccentrocity at reference epoch */
    REAL8 *meanPerAno,                     /**< [out] mean anomaly of periastron */
    REAL8 *deltaT,                         /**< [out] sampling interval (s) */
    REAL8 *f_min,                          /**< [out] starting GW frequency (Hz) */
    REAL8 *f_ref,                          /**< [out] reference frequency (Hz) */
    LALDict *params                        /**< Input lal dictionary with waveform (and optional) parameters **/
)
{
    *m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    *m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    *S1x = XLALSimInspiralWaveformParamsLookupSpin1x(params);
    *S1y = XLALSimInspiralWaveformParamsLookupSpin1y(params);
    *S1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    *S2x = XLALSimInspiralWaveformParamsLookupSpin2x(params);
    *S2y = XLALSimInspiralWaveformParamsLookupSpin2y(params);
    *S2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);
    *distance = XLALSimInspiralWaveformParamsLookupDistance(params);
    *inclination = XLALSimInspiralWaveformParamsLookupInclination(params);
    *phiRef = XLALSimInspiralWaveformParamsLookupRefPhase(params);
    *longAscNodes = XLALSimInspiralWaveformParamsLookupLongAscNodes(params);
    *eccentricity = XLALSimInspiralWaveformParamsLookupEccentricity(params);
    *meanPerAno = XLALSimInspiralWaveformParamsLookupMeanPerAno(params);
    *deltaT = XLALSimInspiralWaveformParamsLookupDeltaT(params);
    *f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    *f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);

    return;
}

/**
 * Insert all the input arguments needed by XLALSimInspiralChooseTDModes() into a laldictionary.
 */
void XLALSimInspiralParseDictionaryToChooseTDModes(
    REAL8 *phiRef,                         /**< [out] reference orbital phase (rad) */
    REAL8 *deltaT,                         /**< [out] sampling interval (s) */
    REAL8 *m1,                             /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                             /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                            /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                            /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                            /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                            /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                            /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                            /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *f_min,                          /**< [out] starting GW frequency (Hz) */
    REAL8 *f_ref,                          /**< [out] reference frequency (Hz) */
    REAL8 *distance,                       /**< [out] distance of source (m) */
    INT4 *lmax,                            /**< [out] generate all modes with l <= lmax */
    LALDict *params                        /**< Input lal dictionary with waveform (and optional) parameters **/
)
{
    *phiRef = XLALSimInspiralWaveformParamsLookupRefPhase(params);
    *deltaT = XLALSimInspiralWaveformParamsLookupDeltaT(params);
    *m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    *m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    *S1x = XLALSimInspiralWaveformParamsLookupSpin1x(params);
    *S1y = XLALSimInspiralWaveformParamsLookupSpin1y(params);
    *S1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    *S2x = XLALSimInspiralWaveformParamsLookupSpin2x(params);
    *S2y = XLALSimInspiralWaveformParamsLookupSpin2y(params);
    *S2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);
    *f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    *f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);
    *distance = XLALSimInspiralWaveformParamsLookupDistance(params);
    *lmax = XLALSimInspiralWaveformParamsLookupLmax(params);

    return;
}

/**
 * Insert all the input arguments needed by XLALSimInspiralChooseFDWaveform() into a laldictionary.
 */
void XLALSimInspiralParseDictionaryToChooseFDWaveform(
    REAL8 *m1,                             /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                             /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                            /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                            /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                            /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                            /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                            /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                            /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *distance,                       /**< [out] distance of source (m) */
    REAL8 *inclination,                    /**< [out] inclination of source (rad) */
    REAL8 *phiRef,                         /**< [out] reference orbital phase (rad) */
    REAL8 *longAscNodes,                   /**< [out] longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 *eccentricity,                   /**< [out] eccentrocity at reference epoch */
    REAL8 *meanPerAno,                     /**< [out] mean anomaly of periastron */
    REAL8 *deltaF,                         /**< [out] frequency interval (Hz) */
    REAL8 *f_min,                          /**< [out] starting GW frequency (Hz) */
    REAL8 *f_max,                          /**< [out] ending GW frequency (Hz) */
    REAL8 *f_ref,                          /**< [out] reference frequency (Hz) */
    LALDict *params                        /**< Input lal dictionary with waveform (and optional) parameters **/
)
{

    *m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    *m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    *S1x = XLALSimInspiralWaveformParamsLookupSpin1x(params);
    *S1y = XLALSimInspiralWaveformParamsLookupSpin1y(params);
    *S1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    *S2x = XLALSimInspiralWaveformParamsLookupSpin2x(params);
    *S2y = XLALSimInspiralWaveformParamsLookupSpin2y(params);
    *S2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);
    *distance = XLALSimInspiralWaveformParamsLookupDistance(params);
    *inclination = XLALSimInspiralWaveformParamsLookupInclination(params);
    *phiRef = XLALSimInspiralWaveformParamsLookupRefPhase(params);
    *longAscNodes = XLALSimInspiralWaveformParamsLookupLongAscNodes(params);
    *eccentricity = XLALSimInspiralWaveformParamsLookupEccentricity(params);
    *meanPerAno = XLALSimInspiralWaveformParamsLookupMeanPerAno(params);
    *deltaF = XLALSimInspiralWaveformParamsLookupDeltaF(params);
    *f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    *f_max = XLALSimInspiralWaveformParamsLookupFMax(params);
    *f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);

    return;
}

/**
 * Insert all the input arguments needed by XLALSimInspiralChooseFDModes() into a laldictionary.
 */
void XLALSimInspiralParseDictionaryToChooseFDModes(
    REAL8 *m1,                                   /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                                   /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                                  /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                                  /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                                  /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                                  /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                                  /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                                  /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *deltaF,                               /**< [out] sampling interval (s) */
    REAL8 *f_min,                                /**< [out] starting GW frequency (Hz) */
    REAL8 *f_max,                                /**< [out] ending GW frequency (Hz) */
    REAL8 *f_ref,                                /**< [out] reference GW frequency (Hz) */
    REAL8 *phiRef,                               /**< [out] reference phase (rad) */
    REAL8 *distance,                             /**< [out] distance of source (m) */
    REAL8 *inclination,                          /**< [out] inclination of source (rad) */
    LALDict *params                              /**< Input lal dictionary with waveform (and optional) parameters **/
)
{
    *m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    *m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    *S1x = XLALSimInspiralWaveformParamsLookupSpin1x(params);
    *S1y = XLALSimInspiralWaveformParamsLookupSpin1y(params);
    *S1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    *S2x = XLALSimInspiralWaveformParamsLookupSpin2x(params);
    *S2y = XLALSimInspiralWaveformParamsLookupSpin2y(params);
    *S2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);
    *deltaF = XLALSimInspiralWaveformParamsLookupDeltaF(params);
    *f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    *f_max = XLALSimInspiralWaveformParamsLookupFMax(params);
    *f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);
    *phiRef = XLALSimInspiralWaveformParamsLookupRefPhase(params);
    *distance = XLALSimInspiralWaveformParamsLookupDistance(params);
    *inclination = XLALSimInspiralWaveformParamsLookupInclination(params);

    return;
}




 
 /** @} */

/**
 * @name General Waveform Switching Generation Routines
 * @{
 */

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the time domain.
 *
 * The parameters passed must be in SI units.
 */
int XLALSimInspiralChooseTDWaveform(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    const REAL8 m1,                             /**< mass of companion 1 (kg) */
    const REAL8 m2,                             /**< mass of companion 2 (kg) */
    const REAL8 S1x,                            /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                            /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                            /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                            /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                            /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                            /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                       /**< distance of source (m) */
    const REAL8 inclination,                    /**< inclination of source (rad) */
    const REAL8 phiRef,                         /**< reference orbital phase (rad) */
    const REAL8 longAscNodes,                   /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 eccentricity,                   /**< eccentrocity at reference epoch */
    const REAL8 UNUSED meanPerAno,              /**< mean anomaly of periastron */
    const REAL8 deltaT,                         /**< sampling interval (s) */
    const REAL8 f_min,                          /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *params,                            /**< LAL dictionary containing accessory parameters */
    const Approximant approximant               /**< post-Newtonian approximant to use for waveform production */
    )
{
	LALSimInspiralGenerator *generator;
    int ret;

    generator = XLALSimInspiralChooseGenerator(approximant, params);
    XLAL_CHECK(generator, XLAL_EFUNC);

    if (params == NULL) {
        params = XLALCreateDict();
    }
    else {
        params = XLALDictDuplicate(params);
    }
    XLAL_CHECK(params, XLAL_EFUNC);

    XLALSimInspiralWaveformParamsInsertMass1(params, m1);
    XLALSimInspiralWaveformParamsInsertMass2(params, m2);
    XLALSimInspiralWaveformParamsInsertSpin1x(params, S1x);
    XLALSimInspiralWaveformParamsInsertSpin1y(params, S1y);
    XLALSimInspiralWaveformParamsInsertSpin1z(params, S1z);
    XLALSimInspiralWaveformParamsInsertSpin2x(params, S2x);
    XLALSimInspiralWaveformParamsInsertSpin2y(params, S2y);
    XLALSimInspiralWaveformParamsInsertSpin2z(params, S2z);
    XLALSimInspiralWaveformParamsInsertDistance(params, distance);
    XLALSimInspiralWaveformParamsInsertInclination(params, inclination);
    XLALSimInspiralWaveformParamsInsertRefPhase(params, phiRef);
    XLALSimInspiralWaveformParamsInsertLongAscNodes(params, longAscNodes);
    XLALSimInspiralWaveformParamsInsertEccentricity(params, eccentricity);
    XLALSimInspiralWaveformParamsInsertMeanPerAno(params, meanPerAno);
    XLALSimInspiralWaveformParamsInsertDeltaT(params, deltaT);
    XLALSimInspiralWaveformParamsInsertF22Start(params, f_min);
    XLALSimInspiralWaveformParamsInsertF22Ref(params, f_ref);

    ret = XLALSimInspiralGenerateTDWaveform(hplus, hcross, params, generator);
    XLALDestroyDict(params);
    XLALDestroySimInspiralGenerator(generator);

    return ret;
}

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the frequency domain.
 */
int XLALSimInspiralChooseFDWaveform(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    const REAL8 m1,                         /**< mass of companion 1 (kg) */
    const REAL8 m2,                         /**< mass of companion 2 (kg) */
    const REAL8 S1x,                        /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                        /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                        /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                        /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                        /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                        /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                   /**< distance of source (m) */
    const REAL8 inclination,                /**< inclination of source (rad) */
    const REAL8 phiRef,                     /**< reference orbital phase (rad) */
    const REAL8 longAscNodes,               /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 eccentricity,               /**< eccentricity at reference epoch */
    const REAL8 UNUSED meanPerAno,          /**< mean anomaly of periastron */
    // frequency sampling parameters, no default value
    const REAL8 deltaF,                     /**< sampling interval (Hz) */
    const REAL8 f_min,                      /**< starting GW frequency (Hz) */
    const REAL8 f_max,                      /**< ending GW frequency (Hz) */
    REAL8 f_ref,                            /**< Reference frequency (Hz) */
    LALDict *params,                     /**< LAL dictionary containing accessory parameters */
    const Approximant approximant           /**< post-Newtonian approximant to use for waveform production */
    )
{
	LALSimInspiralGenerator *generator;
    int ret;

    generator = XLALSimInspiralChooseGenerator(approximant, params);
    XLAL_CHECK(generator, XLAL_EFUNC);

    if (params == NULL) {
        params = XLALCreateDict();
    }
    else {
        params = XLALDictDuplicate(params);
    }
    XLAL_CHECK(params, XLAL_EFUNC);

    /* Avoid duplication of arguments when called from XLALSimInspiralChoose(Generate)TDWaveform  */
    const char *remove_keys[12] = {"total_mass", "chirp_mass", "mass_difference", "reduced_mass", "mass_ratio", "sym_mass_ratio", \
       "spin1_norm", "spin1_tilt", "spin1_phi", "spin2_norm", "spin2_tilt", "spin2_phi"};
    for(size_t j = 0; j < sizeof(remove_keys)/sizeof(*remove_keys); ++j){
      if (XLALDictContains(params, remove_keys[j]))
        XLALDictPop(params, remove_keys[j]);
    }

    XLALSimInspiralWaveformParamsInsertMass1(params, m1);
    XLALSimInspiralWaveformParamsInsertMass2(params, m2);
    XLALSimInspiralWaveformParamsInsertSpin1x(params, S1x);
    XLALSimInspiralWaveformParamsInsertSpin1y(params, S1y);
    XLALSimInspiralWaveformParamsInsertSpin1z(params, S1z);
    XLALSimInspiralWaveformParamsInsertSpin2x(params, S2x);
    XLALSimInspiralWaveformParamsInsertSpin2y(params, S2y);
    XLALSimInspiralWaveformParamsInsertSpin2z(params, S2z);
    XLALSimInspiralWaveformParamsInsertDistance(params, distance);
    XLALSimInspiralWaveformParamsInsertInclination(params, inclination);
    XLALSimInspiralWaveformParamsInsertRefPhase(params, phiRef);
    XLALSimInspiralWaveformParamsInsertLongAscNodes(params, longAscNodes);
    XLALSimInspiralWaveformParamsInsertEccentricity(params, eccentricity);
    XLALSimInspiralWaveformParamsInsertMeanPerAno(params, meanPerAno);
    XLALSimInspiralWaveformParamsInsertDeltaF(params, deltaF);
    XLALSimInspiralWaveformParamsInsertF22Start(params, f_min);
    XLALSimInspiralWaveformParamsInsertFMax(params, f_max);
    XLALSimInspiralWaveformParamsInsertF22Ref(params, f_ref);

    ret = XLALSimInspiralGenerateFDWaveform(hptilde, hctilde, params, generator);
    XLALDestroyDict(params);
    XLALDestroySimInspiralGenerator(generator);

    return ret;
}

/**
 * @brief Generates an time domain inspiral waveform using the specified approximant; the
 * resulting waveform is appropriately conditioned, suitable for injection into data,
 * and decomposed into the (2, \f$\pm\f$ 2), spin -2 weighted spherical harmonic modes.
 * NOTE: This is an algebraic decomposition, and will only be correct for approximants
 * which use only the dominant 2, \f$\pm\f$ 2 mode.
 *
 * For spinning waveforms, all known spin effects up to given PN order are included
 *
 * This routine can generate FD approximants and transform them into the time domain.
 * Waveforms are generated beginning at a slightly lower starting frequency and tapers
 * are put in this early region so that the waveform smoothly turns on.  Artifacts at
 * the very end of the waveform are also tapered.  The resulting waveform is high-pass
 * filtered at frequency f_min so that it should have little content at lower frequencies.
 *
 * This routine used to have one additional parameter relative to XLALSimInspiralChooseTDWaveform:
 * the redshift, z, of the waveform, which is now stuffed into the LALDict structure.
 * This should be set to zero (default value) for sources in the nearby universe (that are nearly at rest relative to the
 * earth).  For sources at cosmological distances, the mass parameters m1 and m2 should
 * be interpreted as the physical (source frame) masses of the bodies and the distance
 * parameter r is the comoving (transverse) distance.  If the calling routine has already
 * applied cosmological "corrections" to m1 and m2 and regards r as a luminosity distance
 * then the redshift factor should again be set to zero.
 *
 * @note The parameters passed must be in SI units.
 */
SphHarmTimeSeries* XLALSimInspiralTDModesFromPolarizations(
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 distance,                             /**< distance of source (m) */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
    REAL8 meanPerAno,                           /**< mean anomaly of periastron */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{

    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) ) {
      XLALPrintError("Non-zero transverse spins were given, but it is not possible to recover modes from H+ and Hx for precessing waveforms.\n");
      XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    REAL8TimeSeries *hplus = NULL, *hcross = NULL;
    COMPLEX16TimeSeries *h22,*h2m2;
    SphHarmTimeSeries *hlm;
    UINT4 j;
    int retval;
    float fac = XLALSpinWeightedSphericalHarmonic(0., 0., -2, 2,2);

    /* Generate waveform via on-axis emission. Assumes only (2,2) and (2,-2) emission */
    retval = XLALSimInspiralTD(&hplus, &hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, 0., phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, LALparams, approximant);
    if (retval < 0)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    /* Step 1: Create COMPLEX16 TimeSeries and populate them */
    h22 = XLALCreateCOMPLEX16TimeSeries("h22", &(hplus)->epoch, 0, deltaT, &lalStrainUnit, (hplus)->data->length);
    h2m2 = XLALCreateCOMPLEX16TimeSeries("h2m2", &(hplus)->epoch, 0, deltaT, &lalStrainUnit, (hplus)->data->length);
    for (j=0; j< (hplus)->data->length; j++) {
      h22->data->data[j] = ((hplus)->data->data[j] - I*((hcross)->data->data[j]))/fac;
      h2m2->data->data[j] = ((hplus)->data->data[j] + I*((hcross)->data->data[j]))/fac;
    }

    /* Step 2: Add them into the data */
    hlm = XLALSphHarmTimeSeriesAddMode(NULL, h22, 2, 2);
    hlm = XLALSphHarmTimeSeriesAddMode(hlm, h2m2, 2, -2);

    /* Step 3: Clean up */
    XLALDestroyREAL8TimeSeries(hplus);
    XLALDestroyREAL8TimeSeries(hcross);
    XLALDestroyCOMPLEX16TimeSeries(h22);
    XLALDestroyCOMPLEX16TimeSeries(h2m2);

    return hlm;
}

/** Helper routines for XLALSimInspiralTD(): performs conditioning of a TD waveform */
int XLALSimInspiralTDFromTD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 distance,                             /**< distance of source (m) */
    REAL8 inclination,                          /**< inclination of source (rad) */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
    REAL8 meanPerAno,                           /**< mean anomaly of periastron */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
)
{
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
    double tchirp, tmerge, textra;
    double fisco, fstart;
    double s;
    int retval;

    if (!XLALSimInspiralImplementedTDApproximants(approximant))
        XLAL_ERROR(XLAL_EINVAL, "Invalid approximant: not a TD approximant");

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    /* apply redshift correction to dimensionful source-frame quantities */
    REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(LALparams);
    if (z != 0.0) {
        m1 *= (1.0 + z);
        m2 *= (1.0 + z);
        distance *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
    }
    /* set redshift to zero so we don't accidentally apply it again later */
    z=0.;
    if (LALparams)
      XLALSimInspiralWaveformParamsInsertRedshift(LALparams,z);

    /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);

    /* upper bound on the final black hole spin */
    s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);

    /* upper bound on the final plunge, merger, and ringdown time */
    tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
    textra = extra_cycles / f_min;

    /* time domain approximant: condition by generating a waveform
     * with a lower starting frequency and apply tapers in the
     * region between that lower frequency and the requested
     * frequency f_min; here compute a new lower frequency */
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp + tmerge + textra, m1, m2);

    /* generate the waveform in the time domain starting at fstart */
    retval = XLALSimInspiralChooseTDWaveform(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, fstart, f_ref, LALparams, approximant);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* condition the time domain waveform by tapering in the extra time
        * at the beginning and high-pass filtering above original f_min */
    XLALSimInspiralTDConditionStage1(*hplus, *hcross, extra_time_fraction * tchirp + textra, original_f_min);

    /* final tapering at the beginning and at the end to remove filter transients */

    /* waveform should terminate at a frequency >= Schwarzschild ISCO
        * so taper one cycle at this frequency at the end; should not make
        * any difference to IMR waveforms */
    fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);

    return 0;
}

/** Helper routines for XLALSimInspiralTD(): performs conditioning of a FD waveform and transforms it to TD */
int XLALSimInspiralTDFromFD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 distance,                             /**< distance of source (m) */
    REAL8 inclination,                          /**< inclination of source (rad) */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
    REAL8 meanPerAno,                           /**< mean anomaly of periastron */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
)
{
    COMPLEX16FrequencySeries *hptilde = NULL;
    COMPLEX16FrequencySeries *hctilde = NULL;
    REAL8FFTPlan *plan;
    size_t chirplen, end, k;
    double tshift;
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
    double f_max = 0.5 / deltaT;
    double tchirp, tmerge, textra;
    double fisco, fstart;
    double s;
    int retval;

    if (!XLALSimInspiralImplementedFDApproximants(approximant))
        XLAL_ERROR(XLAL_EINVAL, "Invalid approximant: not a FD approximant");

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    /* apply redshift correction to dimensionful source-frame quantities */
    REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(LALparams);
    if (z != 0.0) {
        m1 *= (1.0 + z);
        m2 *= (1.0 + z);
        distance *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
    }
    /* set redshift to zero so we don't accidentally apply it again later */
    z=0.;
    if (LALparams)
      XLALSimInspiralWaveformParamsInsertRedshift(LALparams,z);

    /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);

    /* upper bound on the final black hole spin */
    s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);

    /* upper bound on the final plunge, merger, and ringdown time */
    tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
    textra = extra_cycles / f_min;

    /* generate the conditioned waveform in the frequency domain */
    /* note: redshift factor has already been applied above */
    /* set deltaF = 0 to get a small enough resolution */
    retval = XLALSimInspiralFD(&hptilde, &hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, 0.0, f_min, f_max, f_ref, LALparams, approximant);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* we want to make sure that this waveform will give something
     * sensible if it is later transformed into the time domain:
     * to avoid the end of the waveform wrapping around to the beginning,
     * we shift waveform backwards in time and compensate for this
     * shift by adjusting the epoch -- note that XLALSimInspiralFD
     * guarantees that there is extra padding to do this */
    tshift = round(textra / deltaT) * deltaT; /* integer number of samples */
    for (k = 0; k < hptilde->data->length; ++k) {
        double complex phasefac = cexp(2.0 * M_PI * I * k * hptilde->deltaF * tshift);
        hptilde->data->data[k] *= phasefac;
        hctilde->data->data[k] *= phasefac;
    }
    XLALGPSAdd(&hptilde->epoch, tshift);
    XLALGPSAdd(&hctilde->epoch, tshift);

    /* transform the waveform into the time domain */
    chirplen = 2 * (hptilde->data->length - 1);
    *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
    *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
    plan = XLALCreateReverseREAL8FFTPlan(chirplen, 0);
    if (!(*hplus) || !(*hcross) || !plan) {
        XLALDestroyCOMPLEX16FrequencySeries(hptilde);
        XLALDestroyCOMPLEX16FrequencySeries(hctilde);
        XLALDestroyREAL8TimeSeries(*hcross);
        XLALDestroyREAL8TimeSeries(*hplus);
        XLALDestroyREAL8FFTPlan(plan);
        XLAL_ERROR(XLAL_EFUNC);
    }
    XLALREAL8FreqTimeFFT(*hplus, hptilde, plan);
    XLALREAL8FreqTimeFFT(*hcross, hctilde, plan);

    /* apply time domain filter at original f_min */
    XLALHighPassREAL8TimeSeries(*hplus, original_f_min, 0.99, 8);
    XLALHighPassREAL8TimeSeries(*hcross, original_f_min, 0.99, 8);

    /* compute how long a chirp we should have */
    /* revised estimate of chirp length from new start frequency */
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);
    tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, S1z, S2z);

    /* total expected chirp length includes merger */
    chirplen = round((tchirp + tmerge) / deltaT);

    /* amount to snip off at the end is tshift */
    end = (*hplus)->data->length - round(tshift / deltaT);

    /* snip off extra time at beginning and at the end */
    XLALResizeREAL8TimeSeries(*hplus, end - chirplen, chirplen);
    XLALResizeREAL8TimeSeries(*hcross, end - chirplen, chirplen);

    /* clean up */
    XLALDestroyREAL8FFTPlan(plan);
    XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(hctilde);

    /* final tapering at the beginning and at the end to remove filter transients */

    /* waveform should terminate at a frequency >= Schwarzschild ISCO
     * so taper one cycle at this frequency at the end; should not make
     * any difference to IMR waveforms */
    fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);

    return 0;
}

/**
 * @brief Generates an time domain inspiral waveform using the specified approximant; the
 * resulting waveform is appropriately conditioned and suitable for injection into data.
 *
 * For spinning waveforms, all known spin effects up to given PN order are included
 *
 * This routine can generate FD approximants and transform them into the time domain.
 * Waveforms are generated beginning at a slightly lower starting frequency and tapers
 * are put in this early region so that the waveform smoothly turns on.  Artifacts at
 * the very end of the waveform are also tapered.  The resulting waveform is high-pass
 * filtered at frequency f_min so that it should have little content at lower frequencies.

 * If calling with precessing time-domain approximants for which the reference frequency
 * is the starting frequency, or if calling with NR_hdf5 approximant, the starting
 * frequency is not altered. Uses XLALSimInspiralGetSpinFreqFromApproximant to determine
 * appropriate behaviour.
 * Similarly, if calling time-domain models for which a starting frequency of
 * zero is allowed (as set in
 * XLALSimInspiralGetAllowZeroMinFreqFromApproximant), the starting frequency
 * is never altered, independent of f_min.
 *
 * This routine used to have one additional parameter relative to XLALSimInspiralChooseTDWaveform:
 * the redshift, z, of the waveform, which is now stuffed into the LALDict structure.
 * This should be set to zero (default value) for sources in the nearby universe (that are nearly at rest relative to the
 * earth).  For sources at cosmological distances, the mass parameters m1 and m2 should
 * be interpreted as the physical (source frame) masses of the bodies and the distance
 * parameter r is the comoving (transverse) distance.  If the calling routine has already
 * applied cosmological "corrections" to m1 and m2 and regards r as a luminosity distance
 * then the redshift factor should again be set to zero.
 *
 * @note The parameters passed must be in SI units.
 */
int XLALSimInspiralTD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 distance,                             /**< distance of source (m) */
    REAL8 inclination,                          /**< inclination of source (rad) */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
    REAL8 meanPerAno,                           /**< mean anomaly of periastron */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    int ret;

    /* set condition flag */
    if (LALparams == NULL) {
        LALparams = XLALCreateDict();
    } else {
        LALparams = XLALDictDuplicate(LALparams);
        if (XLALDictContains(LALparams, "condition"))
            XLALDictRemove(LALparams, "condition");
    }
    XLALDictInsertINT4Value(LALparams, "condition", 2);

    ret = XLALSimInspiralChooseTDWaveform(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, LALparams, approximant);
    XLALDestroyDict(LALparams);

    return ret;
}

/**
 * @brief Generates a frequency domain inspiral waveform using the specified approximant; the
 * resulting waveform is appropriately conditioned and suitable for injection into data.
 *
 * For spinning waveforms, all known spin effects up to given PN order are included.
 *
 * This routine can generate TD approximants and transform them into the frequency domain.
 * Waveforms are generated beginning at a slightly lower starting frequency and tapers
 * are put in this early region so that the waveform smoothly turns on.
 *
 * If an FD approximant is used, this routine applies tapers in the frequency domain
 * between the slightly-lower frequency and the requested f_min.  Also, the phase of the
 * waveform is adjusted to introduce a time shift.  This time shift should allow the
 * resulting waveform to be Fourier transformed into the time domain without wrapping
 * the end of the waveform to the beginning.
 *
 * This routine assumes that f_max is the Nyquist frequency of a corresponding time-domain
 * waveform, so that deltaT = 0.5 / f_max.  If deltaF is set to 0 then this routine computes
 * a deltaF that is small enough to represent the Fourier transform of a time-domain waveform
 * while ensuring the length of the time-domain signal if a power of 2.
 * If deltaF is specified but f_max / deltaF is not a power of 2, and the waveform approximant
 * is a time-domain approximant, then f_max is increased so that f_max / deltaF is the next
 * power of 2.  (If the user wishes to discard the extra high frequency content, this must
 * be done separately.)
 *
 * The user should take care to ensure that deltaF is sufficiently small to contain the full
 * signal (time series duration = 1 / deltaF). If the provided deltaF is too large the signal
 * will be abruptly truncated for time-domain waveform generators. For frequency-domain
 * generators the signal will be aliased in the frequency domain.
 *
 * Similarly, if the provided f_max is less than the ringdown frequency the underlying waveform
 * generator may raise an error. If not, the frequency domain signal will be aliased in the
 * frequency domain.
 *
 * Some waveform approximants have built in checks for the maximum frequency and signal length
 *
 * This routine used to have one additional parameter relative to XLALSimInspiralChooseTDWaveform:
 * the redshift, z, of the waveform, which is now stuffed into the LALDict.
 * This should be set to zero (default value) for sources in the nearby universe (that are nearly at rest relative to the
 * earth).  For sources at cosmological distances, the mass parameters m1 and m2 should
 * be interpreted as the physical (source frame) masses of the bodies and the distance
 * parameter r is the comoving (transverse) distance.  If the calling routine has already
 * applied cosmological "corrections" to m1 and m2 and regards r as a luminosity distance
 * then the redshift factor should again be set to zero.
 *
 *
 * @note The parameters passed must be in SI units.
 */
int XLALSimInspiralFD(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    REAL8 m1,                               /**< mass of companion 1 (kg) */
    REAL8 m2,                               /**< mass of companion 2 (kg) */
    REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
    REAL8 distance,                         /**< distance of source (m) */
    REAL8 inclination,                      /**< inclination of source (rad) */
    REAL8 phiRef,                           /**< reference orbital phase (rad) */
    REAL8 longAscNodes,                     /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 eccentricity,                     /**< eccentricity at reference epoch */
    REAL8 meanPerAno,                       /**< mean anomaly of periastron */
    REAL8 deltaF,                           /**< sampling interval (Hz) */
    REAL8 f_min,                            /**< starting GW frequency (Hz) */
    REAL8 f_max,                            /**< ending GW frequency (Hz) */
    REAL8 f_ref,                            /**< Reference frequency (Hz) */
    LALDict *LALparams,                     /**< LAL dictionary containing accessory parameters */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
    )
{
    int ret;

    /* set condition flag */
    if (LALparams == NULL) {
        LALparams = XLALCreateDict();
    } else {
        LALparams = XLALDictDuplicate(LALparams);
        if (XLALDictContains(LALparams, "condition"))
            XLALDictRemove(LALparams, "condition");
    }
    XLALDictInsertINT4Value(LALparams, "condition", 2);

    ret = XLALSimInspiralChooseFDWaveform(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, f_min, f_max, f_ref, LALparams, approximant);
    XLALDestroyDict(LALparams);

    return ret;
}

/**
 * @deprecated Use XLALSimInspiralChooseTDWaveform() instead
 *
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 *
 * The parameters passed must be in SI units.
 */
int XLALSimInspiralChooseWaveform(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    const REAL8 m1,                                   /**< mass of companion 1 (kg) */
    const REAL8 m2,                                   /**< mass of companion 2 (kg) */
    const REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                             /**< distance of source (m) */
    const REAL8 inclination,                          /**< inclination of source (rad) */
    const REAL8 phiRef,                               /**< reference orbital phase (rad) */
    const REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
    const REAL8 meanPerAno,                           /**< mean anomaly of periastron */
    // frequency sampling parameters, no default value
    const REAL8 deltaT,                               /**< sampling interval (s) */
    const REAL8 f_min,                                /**< starting GW frequency (Hz) */
    const REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    LALDict *LALpars,                                 /**< LAL dictionary containing accessory parameters */
    const Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralChooseTDWaveform");

    return XLALSimInspiralChooseTDWaveform(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z,
                       distance, inclination, phiRef, longAscNodes,
                       eccentricity, meanPerAno, deltaT, f_min, f_ref,
                       LALpars, approximant);
}

/** @} */

/**
 * @name General Waveform Switching Mode Generation Routines
 * @{
 */

/**
 * Interface to compute a set of -2 spin-weighted spherical harmonic modes
 * for a binary inspiral for a given waveform approximant.
 * PN Approximants (TaylorT1 - T4), EOBNRv2 (EOBNRv2HM), NRSur7dq2, NRSur7dq4
 * NRHybSur3dq8 and spin-precessing SpintaylorT1, T5, T4 are implemented.
 *
 * The EOBNRv2 model returns the (2,2), (2,1), (3,3), (4,4), and (5,5) modes.
 * Note that the inclination parameter is not passed to create hlm modes,
 * hence to recover the correct h+,x one has to combine the hlm modes with
 * Euler angles alpha=0, iota=inclination, psi=0,Pi/2 (according to the approximat) i.e.
 * (h+ + I hx) (psi,iota,alpha)= e^(-2Ialpha) Sum_{l,m} Y_lm(-iota,-psi) h_lm,
 * or equivalently rotate h_lm -> h'_lm=DWigner(psi,iota,alpha) h_lm
 * and then obtain
 * (h+ + I hx) = Sum_{l,m} Y_lm(0,0) h'_lm,
 */
SphHarmTimeSeries *XLALSimInspiralChooseTDModes(
    UNUSED REAL8 phiRef,                        /**< reference orbital phase (rad). This variable is not used and only kept here for backwards compatibility */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 distance,                             /**< distance of source (m) */
    LALDict *params,                            /**< LAL dictionary containing accessory parameters */
    int lmax,                                   /**< generate all modes with l <= lmax */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
	LALSimInspiralGenerator *generator;
    SphHarmTimeSeries *hlms = NULL;

    generator = XLALSimInspiralChooseGenerator(approximant, params);
	if (!(generator))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    if (params == NULL) {
        params = XLALCreateDict();
    }
    else {
        params = XLALDictDuplicate(params);
    }
	if (!(params))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    XLALSimInspiralWaveformParamsInsertRefPhase(params, phiRef);
    XLALSimInspiralWaveformParamsInsertDeltaT(params, deltaT);
    XLALSimInspiralWaveformParamsInsertMass1(params, m1);
    XLALSimInspiralWaveformParamsInsertMass2(params, m2);
    XLALSimInspiralWaveformParamsInsertSpin1x(params, S1x);
    XLALSimInspiralWaveformParamsInsertSpin1y(params, S1y);
    XLALSimInspiralWaveformParamsInsertSpin1z(params, S1z);
    XLALSimInspiralWaveformParamsInsertSpin2x(params, S2x);
    XLALSimInspiralWaveformParamsInsertSpin2y(params, S2y);
    XLALSimInspiralWaveformParamsInsertSpin2z(params, S2z);
    XLALSimInspiralWaveformParamsInsertF22Start(params, f_min);
    XLALSimInspiralWaveformParamsInsertF22Ref(params, f_ref);
    XLALSimInspiralWaveformParamsInsertDistance(params, distance);
    XLALSimInspiralWaveformParamsInsertLmax(params, lmax);

    XLALSimInspiralGenerateTDModes(&hlms, params, generator);
    XLALDestroyDict(params);
    XLALDestroySimInspiralGenerator(generator);

    return hlms;
}

/**
 * @brief Interface to compute a set of -2 spin-weighted spherical harmonic modes
 * for a binary merger for a given waveform approximant in the Fourier domain.
 * Non-precessing models IMRPhenomXHM, SEOBNRv4HM_ROM, SEOBNRv5(HM)_ROM and IMRPhenomHM and the
 * precessing IMRPhenomXPHM are implemented.
 * By default, all the modes available in the model are returned, although the list
 * can be specified through the ModeArray option in the LAL dictionary.
 * @details
 * In the Fourier domain the modes span over the whole frequency regime (positive and negative frequencies).
 * However, in the aligned spin case, the modes have support only in one half of the frequency regime.
 * The LAL conventions establish that the negative modes (m<0) have support for positive frequencies while
 * the positive modes (m>0) have support for negative frequencies (this is based in the right hand rule and
 * Fourier transform definition in LAL).
 * Due to the equatorial symmetry of non-precessng systems, there exist a relation between them: h_{lm}(f) = (-1)^l h*_{l-m}(-f).
 * In the precessing case this symmetry is broken and all the modes have support for both positive and negative frequencies.
 *
 * For this reason, the ouput SphHarmFrequencySeries object will return the modes in the whole frequency range.
 * The frequencies of this object will be sorted as: -f_max,...,-f_min,...0,...,f_min,...,f_max.
 * The values of the waveform will be sorted correspondingly. Consequently, in the aligned spin case, half of the frequency spectrum consists of zeros.
 *
 * It is relevant to mention why the arguments inclination and phiRef are needed for computing the h_lm.
 * For AS models the argument inclination is irrelevant and will not be use since it only enters in the Ylm. However, for the precessing model,
 * since the modes are returned in the J-frame, we need the inclination argument to carry out the Euler transformation from the co-precessing L-frame
 * to the inertial J-frame. Regarding the argument phiRef, this affects the output of the precessing model due to the same reason as before, while for the AS models 
 * it would not affect the output of SEOBNRv4HM_ROM, SEOBNRv5(HM)_ROM but will change the output of IMRPhenomHM and IMRPhenomXHM (this is due to the internals workings of the models).
 * If one wants to built the polarizations from the individual modes of ChooseFDModes must be aware of this behaviour. 
 * Ideally, one would call ChooseFDModes with phiRef=0 to obtain the h_lms, then build the Fourier domain polarizations as
 *
 * h_+ (f) = 1/2 sum_{l=2} sum_{m=-l}^{m=l}  (  h_lm(f) * Y_lm(theta, vphi)  +  h*_lm(-f) * Y*_lm(theta, vphi)  )
 * h_x (f) = i/2 sum_{l=2} sum_{m=-l}^{m=l}  (  h_lm(f) * Y_lm(theta, vphi)  -  h*_lm(-f) * Y*_lm(theta, vphi)  )
 *
 * where theta is the inclination and vphi is pi/2 - phiRef.
 *
 * If one does this, one will find generally a very close result to ChooseFDWaveform with very small mismatches (~10^-9 for IMRPhenomXHM),
 * and this is what is found if one uses the XLALSimInspiralPolarizationsFromSphHarmFrequencySeries function, however this is not close to machine precision.
 * The reason is that IMRPhenomHM and IMRPhenomXHM use internally the phiRef argument to compute the h_lms because at that time phiRef was considered to be also a reference
 * phase for the h_lms and not just the argument for the azimuthal angle in the Y_lms.
 * To take this into account we provide also the function XLALSimInspiralPolarizationsFromChooseFDModes which build the polarizations in the proper way for each
 * model returning a result close to machine precision with ChooseFDWaveform.
 *
 * For the precessing model IMRPhenomXPHM, since the h_lms are returned in the J-frame one must build the polarizations using theta = theta_JN and vphi = 0.
 * The parameter theta_JN is computed internally when using XLALSimInspiralPolarizationsFromChooseFDModes and again here the result is close to machine precision to ChooseFDWaveform.
 * However, one would have to compute theta_JN personally when using XLALSimInspiralPolarizationsFromSphHarmFrequencySeries.
 * For IMRPhenomXPHM this parameter can be computed using XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame. Eventhough one would still have to correct with the polarization angle.
 *
 * By default all the modes available in the model will be returned, both positive and negative modes.
 * The mode content of AS models can be adjusted through the ModeArray option in the LAL dictionary argument, and it accepts any set of modes e.g. (2,2),(2,-1),(3,3),...
 * For IMRPhenomXPHM, we can specify both the modes in the co-precessing L-frame, which are used to do the twisting-up, and
 * the ouput modes in the inertial J-frame. The set of modes in the L-frame are specified with the standard ModeArray option
 * while the set of modes in the J-frame are specified in a new option called ModeArrayJframe. Notice that in IMRPhenomXPHM, ModeArray does not distinguish
 * between positive or negative modes and it always twists-up both positive and negative modes, i.e. the sets {(2,2)}, {(2,-2)} or {(2,2),(2,-2)} would return the same result.
 * For the modes in ModeArrayJframe, we can specify both positive or negative modes for example {(2,2),(2,-2),(2,-1),(3,3),...}.
 */
SphHarmFrequencySeries *XLALSimInspiralChooseFDModes(
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 deltaF,                               /**< sampling interval (s) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_max,                                /**< ending GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 phiRef,								/**< reference phase (rad) */
    REAL8 distance,                             /**< distance of source (m) */
    REAL8 inclination,                          /**< inclination of source (rad) */
	LALDict *params,                            /**< LAL dictionary containing accessory parameters (optional mode array) */
    Approximant approximant                     /**< approximant to use for waveform production */
    )
{

	LALSimInspiralGenerator *generator;
    SphHarmFrequencySeries *hlms = NULL;

    generator = XLALSimInspiralChooseGenerator(approximant, params);
	if (!(generator))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    if (params == NULL) {
        params = XLALCreateDict();
    }
    else {
        params = XLALDictDuplicate(params);
    }
	if (!(params))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    XLALSimInspiralWaveformParamsInsertMass1(params, m1);
    XLALSimInspiralWaveformParamsInsertMass2(params, m2);
    XLALSimInspiralWaveformParamsInsertSpin1x(params, S1x);
    XLALSimInspiralWaveformParamsInsertSpin1y(params, S1y);
    XLALSimInspiralWaveformParamsInsertSpin1z(params, S1z);
    XLALSimInspiralWaveformParamsInsertSpin2x(params, S2x);
    XLALSimInspiralWaveformParamsInsertSpin2y(params, S2y);
    XLALSimInspiralWaveformParamsInsertSpin2z(params, S2z);
    XLALSimInspiralWaveformParamsInsertDeltaF(params, deltaF);
    XLALSimInspiralWaveformParamsInsertF22Start(params, f_min);
    XLALSimInspiralWaveformParamsInsertFMax(params, f_max);
    XLALSimInspiralWaveformParamsInsertF22Ref(params, f_ref);
    XLALSimInspiralWaveformParamsInsertRefPhase(params, phiRef);
    XLALSimInspiralWaveformParamsInsertDistance(params, distance);
    XLALSimInspiralWaveformParamsInsertInclination(params, inclination);


    XLALSimInspiralGenerateFDModes(&hlms, params, generator);
    XLALDestroyDict(params);
    XLALDestroySimInspiralGenerator(generator);

	return hlms;
}

/**
 * @brief Interface to compute a conditioned set of -2 spin-weighted spherical
 * harmonic modes for a binary inspiral
 * @details
 * This routine is a wrapper for XLALSimInspiralChooseTDModes which applies
 * waveform conditioning to the modes generated by that routine.  The conditioning
 * algorithm is analogous to that performed in XLALSimInspiralTD.  Note that
 * the modes are high-pass filtered at frequency f_min, which is specified for
 * the m=2 mode, which means that the low-frequency part of the m=1 mode is
 * removed by the filtering.  The phasing is computed with any of the TaylorT1,
 * T2, T3, T4 methods.  It can also return the (2,2), (2,1), (3,3), (4,4),
 * (5,5) modes of the EOBNRv2 model. Note that EOBNRv2 will ignore ampO,
 * phaseO, lmax and f_ref arguments.
 * @param deltaT      Sampling interval (s)
 * @param m1          Mass of companion 1 (kg)
 * @param m2          Mass of companion 2 (kg)
 * @param f_min       Starting GW frequency (Hz)
 * @param f_ref       Reference GW frequency (Hz)
 * @param r           Distance of source (m)
 * @param LALpars     LAL dictionary containing accesory parameters
 * @param lmax        Generate all modes with l <= lmax
 * @param approximant Post-Newtonian approximant to use for waveform production
 * @return Linked list of SphHarmTimeSeries modes.
 */
SphHarmTimeSeries *XLALSimInspiralModesTD(REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, LALDict *LALpars, int lmax, Approximant approximant)
{
    const size_t min_taper_samples = 4;
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
    //double tchirp, tmerge, textra;
    //double fisco, fstart;
    double tchirp, textra;
    double fisco;
    //double s;
    size_t length, nzeros, ntaper;
    size_t j;
    SphHarmTimeSeries *modes, *hlm;

    /* if the requested low frequency is below the lowest Kerr ISCO frequency
     * then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, 0.0, 0.0);

    /* upper bound on the final black hole spin */
    //s = XLALSimInspiralFinalBlackHoleSpinBound(0.0, 0.0);

    /* upper bound on the final plunge, merger, and ringdown time */
    //tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
    textra = extra_cycles / f_min;

    /* condition by generating a waveform with a lower starting frequency and
     * apply tapers in the region between that lower frequency and the
     * requested frequency f_min; here compute a new lower frequency */
    // fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp + tmerge + textra, m1, m2);

    XLALPrintWarning("XLAL Warning - XLALSimInspiralModesTD does not yet implement spins - passing zeros\n");
    modes = XLALSimInspiralChooseTDModes(0.,deltaT, m1, m2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, f_min, f_ref, r, LALpars, lmax, approximant);
    if (!modes)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    /* Note: fstart and f_min are frequencies for a m=2 mode, while the
     * frequencies of the mth mode are f_m = m * f / 2; this means that if we
     * generate a waveform that starts at fstart in the m=2 mode then the m=1
     * mode starts at fstart/2 while the m=3 mode starts at 3*fstart/2, and so
     * on.  However, the time it takes the m=2 mode to go from fstart to f_min
     * is the same as the time that a m=1 mode to go from fstart/2 to f_min/2,
     * etc., so we apply tapers over this time.  The upshot is that the
     * resulting modes will be valid for frequencies above m * f_min / 2.
     * evolve from fstart to f_min is the same as the time to evolve from
     * 2*f_start/m to 2*f_min/m for the mth mode, so taper all modes over
     * taper this duration. */
    length = nzeros = modes->mode->data->length;
    for (hlm = modes; hlm; hlm = hlm->next) {
        /* some waveform generators zero-pad the end of the waveform, and we
         * want to remove this; however, we want all modes to have the same
         * length timeseries, so find the minimum number of zeros to excise */
        if (nzeros) {
            j = 0;
            while (hlm->mode->data->data[hlm->mode->data->length - j - 1] == 0.0)
                ++j;
            if (j < nzeros)
                nzeros = j;
        }

        /* here is where we taper the beginning of the waveform below f_min */
        ntaper = round((extra_time_fraction * tchirp + textra) / deltaT);
        for (j = 0; j < ntaper; ++j)
             hlm->mode->data->data[j] *= 0.5 - 0.5 * cos(j * LAL_PI / ntaper);

        /* now high-pass filter the data at the original f_min value so that
         * the modes have negligible frequency content below that frequency;
         * note: this will cut off the low-frequency content of the m=1 mode */
        XLALHighPassCOMPLEX16TimeSeries(hlm->mode, original_f_min, 0.99, 8);
    }

    /* new length after clipping zeros from end */
    length -= nzeros;
    if (nzeros)
        XLALResizeSphHarmTimeSeries(modes, 0, length);

    /* stage 2 conditioning: final tapering at beginning and end */
    /* final tapering at the beginning and at the end */
    /* if this waveform is shorter than 2*min_taper_samples, do nothing */
    if (length < 2 * min_taper_samples) {
        XLAL_PRINT_WARNING("waveform is too shorter than %zu samples: no final tapering applied", 2 * min_taper_samples);
        return modes;
    }

    /* waveform should terminate at a frequency >= Schwarzschild ISCO
     * so taper one cycle at this frequency at the end; should not make
     * any difference to IMR waveforms */
    fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    ntaper = round(1.0 / (fisco * deltaT));
    if (ntaper < min_taper_samples)
        ntaper = min_taper_samples;
    for (hlm = modes; hlm; hlm = hlm->next)
        for (j = 1; j < ntaper; ++j)
             hlm->mode->data->data[length - j] *= 0.5 - 0.5 * cos(j * LAL_PI / ntaper);

    /* there could be a filter transient at the beginning too: we should have
     * some safety there owing to the fact that we are starting at a lower
     * frequency than is needed, so taper off one cycle at the low frequency */
    ntaper = round(1.0 / (f_min * deltaT));
    if (ntaper < min_taper_samples)
        ntaper = min_taper_samples;
    for (hlm = modes; hlm; hlm = hlm->next)
        for (j = 1; j < ntaper; ++j)
             hlm->mode->data->data[j] *= 0.5 - 0.5 * cos(j * LAL_PI / ntaper);

    return modes;
}

/**
 * Interface to compute a single -2 spin-weighted spherical harmonic mode
 * for a binary inspiral of any available amplitude and phase PN order.
 * The phasing is computed with any of the TaylorT1, T2, T3, T4 methods.
 */
COMPLEX16TimeSeries *XLALSimInspiralChooseTDMode(
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    int l,                                      /**< l index of mode */
    int m,                                      /**< m index of mode */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 v0 = 1.;
    COMPLEX16TimeSeries *hlm;
    SphHarmTimeSeries *ts;

    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    if( nonGRparams )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if( deltaT < 1./16385. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    switch (approximant)
    {
        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT1PNMode(v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT2PNMode(v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT3PNMode(v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT4PNMode(v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case EOBNRv2:
        case EOBNRv2HM:
            ts = XLALSimIMREOBNRv2Modes(deltaT, m1, m2, f_min, r);
            hlm = XLALSphHarmTimeSeriesGetMode(ts, l, m);
            break;

        default:
            XLALPrintError("Cannot generate modes for this approximant\n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlm;
}

/** @} */

/**
 * @name Routines for Generating Inspiral Waveforms from Orbital Data
 * @{
 */

/**
 * Given time series for a binary's orbital dynamical variables,
 * construct the waveform polarizations h+ and hx directly.
 * NB: Valid only for non-precessing binaries!
 *
 * Implements Equations (8.8) - (8.10) of:
 * Luc Blanchet, Guillaume Faye, Bala R. Iyer and Siddhartha Sinha,
 * \"The third post-Newtonian gravitational wave polarisations
 * and associated spherical harmonic modes for inspiralling compact binaries
 * in quasi-circular orbits\", Class. Quant. Grav. 25 165003 (2008);
 * arXiv:0802.1249.
 * NB: Be sure to check arXiv:0802.1249v3, which corrects a typo!
 *
 * Note however, that we do not include the constant \"memory\" terms
 */
int XLALSimInspiralPNPolarizationWaveforms(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
    REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
    REAL8TimeSeries *V,       /**< post-Newtonian (PN) parameter */
    REAL8TimeSeries *Phi,     /**< orbital phase */
    REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
    REAL8 m1,                 /**< mass of companion 1 (kg) */
    REAL8 m2,                 /**< mass of companion 2 (kg) */
    REAL8 r,                  /**< distance of source (m) */
    REAL8 i,                  /**< inclination of source (rad) */
    int ampO                  /**< twice PN order of the amplitude */
    )
{
    REAL8 M, eta, eta2, eta3, dm, dist, ampfac, phi, phiShift, v, v2, v3;
    REAL8 hp0, hp05, hp1, hp15, hp2, hp25, hp3;
    REAL8 hc0, hc05, hc1, hc15, hc2, hc25, hc3;
    REAL8 ci, si, ci2, ci4, ci6, ci8, si2, si3, si4, si5, si6;
    INT4 idx, len;

    /* Sanity check input time series */
    LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length
            * sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length
            * sizeof(*(*hcross)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    eta2 = eta*eta;	eta3 = eta2*eta;
    dm = (m1 - m2) / M; // frac. mass difference - \delta m/m in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;

    /*
     * cosines and sines of inclination between
     * line of sight (N) and binary orbital angular momentum (L_N)
     */
    ci = cos(i);  	si = sin(i);
    ci2 = ci*ci;  ci4 = ci2*ci2;  ci6 = ci2*ci4;  ci8 = ci6*ci2;
    si2 = si*si;  si3 = si2*si;  si4 = si2*si2;  si5 = si*si4; si6 = si4*si2;

    /* loop over time steps and compute polarizations h+ and hx */

    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi = Phi->data->data[idx]; 	v = V->data->data[idx];
        v2 = v * v; 	v3 = v * v2;

        /*
         * As explained in Blanchet et al, a phase shift can be applied
         * to make log terms vanish which would appear in the amplitude
         * at 1.5PN and 2.5PN orders. This shift is given in Eq. (8.8)
         * We apply the shift only for the PN orders which need it.
         */
        if( (ampO == -1) || ampO >= 5 )
            phiShift = 3.*v3*(1. - v2*eta/2.)*log( v2 / v0 / v0  );
        else if( ampO >= 3 )
            phiShift = 3.*v3*log( v2 / v0 / v0 );
        else
            phiShift = 0.;

        phi = phi - phiShift;

        /*
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hp0 = hp05 = hp1 = hp15 = hp2 = hp25 = hp3 = 0.;
        hc0 = hc05 = hc1 = hc15 = hc2 = hc25 = hc3 = 0.;

        switch( ampO )
        {
            /* case LAL_PNORDER_THREE_POINT_FIVE: */
            case 7:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            case -1: // Highest known PN order - move if higher terms added!
            /* case LAL_PNORDER_THREE: */
            case 6:
        /* The reference had a typo in the 3PN terms and needed an errata.
         * These should match arXiv:0802.1249v3, which has the fix. */
                hp3 = LAL_PI*dm*si*cos(phi)*(19./64. + ci2*5./16. - ci4/192.
                        + eta*(-19./96. + ci2*3./16. + ci4/96.)) + cos(2.*phi)
                        * (-465497./11025. + (LAL_GAMMA*856./105.
                        - 2.*LAL_PI*LAL_PI/3. + log(16.*v2)*428./105.)
                        * (1. + ci2) - ci2*3561541./88200. - ci4*943./720.
                        + ci6*169./720. - ci8/360. + eta*(2209./360.
                        - LAL_PI*LAL_PI*41./96.*(1. + ci2) + ci2*2039./180.
                        + ci4*3311./720. - ci6*853./720. + ci8*7./360.)
                        + eta2*(12871./540. - ci2*1583./60. - ci4*145./108.
                        + ci6*56./45. - ci8*7./180.) + eta3*(-3277./810.
                        + ci2*19661./3240. - ci4*281./144. - ci6*73./720.
                        + ci8*7./360.)) + LAL_PI*dm*si*cos(3.*phi)*(-1971./128.
                        - ci2*135./16. + ci4*243./128. + eta*(567./64.
                        - ci2*81./16. - ci4*243./64.)) + si2*cos(4.*phi)
                        * (-2189./210. + ci2*1123./210. + ci4*56./9.
                        - ci6*16./45. + eta*(6271./90. - ci2*1969./90.
                        - ci4*1432./45. + ci6*112./45.) + eta2*(-3007./27.
                        + ci2*3493./135. + ci4*1568./45. - ci6*224./45.)
                        + eta3*(161./6. - ci2*1921./90. - ci4*184./45.
                        + ci6*112./45.)) + dm*cos(5.*phi)*(LAL_PI*3125./384.
                        * si3*(1. + ci2)*(1. - 2.*eta)) + si4*cos(6.*phi)
                        * (1377./80. + ci2*891./80. - ci4*729./280.
                        + eta*(-7857./80. - ci2*891./16. + ci4*729./40.)
                        + eta2*(567./4. + ci2*567./10. - ci4*729./20.)
                        + eta3*(-729./16. - ci2*243./80. + ci4*729./40.))
                        + cos(8.*phi)*(-1024./315.*si6*(1. + ci2)*(1. - 7.*eta
                        + 14.*eta2 - 7.*eta3)) + dm*si*sin(phi)*(-2159./40320.
                        - log(2.)*19./32. + (-95./224. - log(2.)*5./8.)*ci2
                        + (181./13440. + log(2.)/96.)*ci4 + eta*(1369./160.
                        + log(2.)*19./48. + (-41./48. - log(2.)*3./8.)*ci2
                        + (-313./480. - log(2.)/48.)*ci4)) + sin(2.*phi)
                        * (-428.*LAL_PI/105.*(1. + ci2)) + dm*si*sin(3.*phi)
                        * (205119./8960. - log(3./2.)*1971./64.
                        + (1917./224. - log(3./2.)*135./8.)*ci2
                        + (-43983./8960. + log(3./2.)*243./64.)*ci4 + eta
                        * (-54869./960. + log(3./2.)*567./32.
                        + (-923./80. - log(3./2.)*81./8.)*ci2
                        + (41851./2880. - log(3./2.)*243./32.)*ci4))
                        + dm*si3*(1. + ci2)*sin(5.*phi)*(-113125./5376.
                        + log(5./2.)*3125./192.
                        + eta*(17639./320. - log(5./2.)*3125./96.));
                hc3 = dm*si*ci*cos(phi)*(11617./20160. + log(2.)*21./16.
                        + (-251./2240. - log(2.)*5./48.)*ci2
                        + eta*(-2419./240. - log(2.)*5./24.
                        + (727./240. + log(2.)*5./24.)*ci2)) + ci*cos(2.*phi)
                        * (LAL_PI*856./105.) + dm*si*ci*cos(3.*phi)
                        * (-36801./896. + log(3./2.)*1809./32.
                        + (65097./4480. - log(3./2.)*405./32.)*ci2
                        + eta*(28445./288. - log(3./2.)*405./16.
                        + (-7137./160. + log(3./2.)*405./16.)*ci2))
                        + dm*si3*ci*cos(5.*phi)*(113125./2688.
                        - log(5./2.)*3125./96. + eta*(-17639./160.
                        + log(5./2.)*3125./48.)) + LAL_PI*dm*si*ci*sin(phi)
                        * (21./32. - ci2*5./96. + eta*(-5./48. + ci2*5./48.))
                        + ci*sin(2.*phi)*(-3620761./44100.
                        + LAL_GAMMA*1712./105. - 4.*LAL_PI*LAL_PI/3.
                        + log(16.*v2)*856./105. - ci2*3413./1260.
                        + ci4*2909./2520. - ci6/45. + eta*(743./90.
                        - 41.*LAL_PI*LAL_PI/48. + ci2*3391./180.
                        - ci4*2287./360. + ci6*7./45.) + eta2*(7919./270.
                        - ci2*5426./135. + ci4*382./45. - ci6*14./45.)
                        + eta3*(-6457./1620. + ci2*1109./180. - ci4*281./120.
                        + ci6*7./45.)) + LAL_PI*dm*si*ci*sin(3.*phi)
                        * (-1809./64. + ci2*405./64. + eta*(405./32.
                        - ci2*405./32.)) + si2*ci*sin(4.*phi)*(-1781./105.
                        + ci2*1208./63. - ci4*64./45. + eta*(5207./45.
                        - ci2*536./5. + ci4*448./45.) + eta2*(-24838./135.
                        + ci2*2224./15. - ci4*896./45.) + eta3*(1703./45.
                        - ci2*1976./45. + ci4*448./45.)) + dm*sin(5.*phi)
                        * (3125.*LAL_PI/192.*si3*ci*(1. - 2.*eta))
                        + si4*ci*sin(6.*phi)*(9153./280. - ci2*243./35.
                        + eta*(-7371./40. + ci2*243./5.) + eta2*(1296./5.
                        - ci2*486./5.) + eta3*(-3159./40. + ci2*243./5.))
                        + sin(8.*phi)*(-2048./315.*si6*ci*(1. - 7.*eta
                        + 14.*eta2 - 7.*eta3));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /* case LAL_PNORDER_TWO_POINT_FIVE: */
            case 5:
                hp25 = cos(phi)*si*dm*(1771./5120. - ci2*1667./5120.
                        + ci4*217./9216. - ci6/9126. + eta*(681./256.
                        + ci2*13./768. - ci4*35./768. + ci6/2304.)
                        + eta2*(-3451./9216. + ci2*673./3072. - ci4*5./9216.
                        - ci6/3072.)) + cos(2.*phi)*LAL_PI*(19./3. + 3.*ci2
                        - ci4*2./3. + eta*(-16./3. + ci2*14./3. + 2.*ci4))
                        + cos(3.*phi)*si*dm*(3537./1024. - ci2*22977./5120.
                        - ci4*15309./5120. + ci6*729./5120.
                        + eta*(-23829./1280. + ci2*5529./1280.
                        + ci4*7749./1280. - ci6*729./1280.)
                        + eta2*(29127./5120. - ci2*27267./5120.
                        - ci4*1647./5120. + ci6*2187./5120.)) + cos(4.*phi)
                        * (-16.*LAL_PI/3.*(1. + ci2)*si2*(1. - 3.*eta))
                        + cos(5.*phi)*si*dm*(-108125./9216. + ci2*40625./9216.
                        + ci4*83125./9216. - ci6*15625./9216.
                        + eta*(8125./256. - ci2*40625./2304. - ci4*48125./2304.
                        + ci6*15625./2304.) + eta2*(-119375./9216.
                        + ci2*40625./3072. + ci4*44375./9216.
                        - ci6*15625./3072.)) + cos(7.*phi)*dm
                        * (117649./46080.*si5*(1. + ci2)*(1. - 4.*eta
                        + 3.*eta2)) + sin(2.*phi)*(-9./5. + ci2*14./5.
                        + ci4*7./5. + eta*(32. + ci2*56./5. - ci4*28./5.))
                        + sin(4.*phi)*si2*(1. + ci2)*(56./5. - 32.*log(2.)/3.
                        + eta*(-1193./30. + 32.*log(2.)));
                /* below would have a constant memory term of si2*ci*eta*6./5. */
                hc25 = cos(2.*phi)*ci*(2. - ci2*22./5. + eta*(-282./5.
                        + ci2*94./5.)) + cos(4.*phi)*ci*si2*(-112./5.
                        + 64.*log(2.)/3. + eta*(1193./15. - 64.*log(2.)))
                        + sin(phi)*si*ci*dm*(-913./7680. + ci2*1891./11520.
                        - ci4*7./4608. + eta*(1165./384. - ci2*235./576.
                        + ci4*7./1152.) + eta2*(-1301./4608. + ci2*301./2304.
                        - ci4*7./1536.)) + sin(2.*phi)*LAL_PI*ci*(34./3.
                        - ci2*8./3. + eta*(-20./3. + 8.*ci2))
                        + sin(3.*phi)*si*ci*dm*(12501./2560. - ci2*12069./1280.
                        + ci4*1701./2560. + eta*(-19581./640. + ci2*7821./320.
                        - ci4*1701./640.) + eta2*(18903./2560.
                        - ci2*11403./1280. + ci4*5103./2560.))
                        + sin(4.*phi)*si2*ci*(-32.*LAL_PI/3.*(1. - 3.*eta))
                        + sin(5.*phi)*si*ci*dm*(-101875./4608. + ci2*6875./256.
                        - ci4*21875./4608. + eta*(66875./1152.
                        - ci2*44375./576. + ci4*21875./1152.)
                        + eta2*(-100625./4608. + ci2*83125./2304.
                        - ci4*21875./1536.)) + sin(7.*phi)*si5*ci*dm
                        * (117649./23040.*(1. - 4.*eta + 3.*eta2));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /* case LAL_PNORDER_TWO: */
            case 4:
                hp2 = cos(phi)*LAL_PI*si*dm*(-5./8. - ci2/8.)
                        + cos(2.*phi)*(11./60. + ci2*33./10. + ci4*29./24.
                        - ci6/24. + eta*(353./36. - 3.*ci2 - ci4*251./72.
                        + ci6*5./24.) + eta2*(-49./12. + ci2*9./2.
                        - ci4*7./24. - ci6*5./24.)) + cos(3.*phi)*LAL_PI*si*dm
                        * (27./8.*(1 + ci2)) + cos(4.*phi)*si2*2./15.*(59.
                        + ci2*35. - ci4*8.
                        - eta*5./3.*(131. + 59.*ci2 + 24.*ci4)
                        + eta2*5.*(21. - 3.*ci2 - 8.*ci4))
                        + cos(6.*phi)*(-81./40.*si4*(1. + ci2)
                        * (1. - 5.*eta + 5.*eta2)) + sin(phi)*si*dm
                        * (11./40. + 5.*log(2)/4. + ci2*(7./40. + log(2)/4.))
                        + sin(3.*phi)*si*dm*((-189./40. + 27./4.*log(3./2.))
                        * (1. + ci2));
                hc2 = cos(phi)*si*ci*dm*(-9./20. - 3./2.*log(2.))
                        + cos(3.*phi)*si*ci*dm*(189./20. - 27./2.*log(3./2.))
                        - sin(phi)*si*ci*dm*3.*LAL_PI/4.
                        + sin(2.*phi)*ci*(17./15. + ci2*113./30. - ci4/4.
                        + eta*(143./9. - ci2*245./18. + ci4*5./4.)
                        + eta2*(-14./3. + ci2*35./6. - ci4*5./4.))
                        + sin(3.*phi)*si*ci*dm*27.*LAL_PI/4.
                        + sin(4.*phi)*ci*si2*4./15.*(55. - 12.*ci2
                        - eta*5./3.*(119. - 36.*ci2)
                        + eta2*5.*(17. - 12.*ci2))
                        + sin(6.*phi)*ci*(-81./20.*si4
                        * (1. - 5.*eta + 5.*eta2));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /* case LAL_PNORDER_ONE_POINT_FIVE: */
            case 3:
                hp15 = cos(phi)*si*dm*(19./64. + ci2*5./16. - ci4/192.
                        + eta*(-49./96. + ci2/8. + ci4/96.))
                        + cos(2.*phi)*(-2.*LAL_PI*(1. + ci2))
                        + cos(3.*phi)*si*dm*(-657./128. - ci2*45./16.
                        + ci4*81./128. + eta*(225./64. - ci2*9./8.
                        - ci4*81./64.)) + cos(5.*phi)*si*dm*(625./384.*si2
                        * (1. + ci2)*(1. - 2.*eta));
                hc15 = sin(phi)*si*ci*dm*(21./32. - ci2*5./96.
                        + eta*(-23./48. + ci2*5./48.))
                        - 4.*LAL_PI*ci*sin(2.*phi) + sin(3.*phi)*si*ci*dm
                        * (-603./64. + ci2*135./64.
                        + eta*(171./32. - ci2*135./32.))
                        + sin(5.*phi)*si*ci*dm*(625./192.*si2*(1. - 2.*eta));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /* case LAL_PNORDER_ONE: */
            case 2:
                hp1 = cos(2.*phi)*(19./6. + 3./2.*ci2 - ci4/3.
                        + eta*(-19./6. + ci2*11./6. + ci4))
                        - cos(4.*phi) * (4./3.*si2*(1. + ci2)*(1. - 3.*eta));
                hc1 = sin(2.*phi)*ci*(17./3. - ci2*4./3.
                        + eta*(-13./3. + 4.*ci2))
                        + sin(4.*phi)*ci*si2*(-8./3.*(1. - 3.*eta));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /*case LAL_PNORDER_HALF:*/
            case 1:
                hp05 = - si*dm*(cos(phi)*(5./8. + ci2/8.)
                        - cos(3.*phi)*(9./8. + 9.*ci2/8.));
                hc05 = si*ci*dm*(-sin(phi)*3./4. + sin(3.*phi)*9./4.);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 0:
                /* below would have a constant memory term of -si2/96.*(17. + ci2) */
                hp0 = -(1. + ci2)*cos(2.*phi);
                hc0 = -2.*ci*sin(2.*phi);
                break;
            /*case LAL_PNORDER_NEWTONIAN:*/
            default:
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d\n",
                        __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
        (*hplus)->data->data[idx] = ampfac * v2 * ( hp0 + v * ( hp05
                + v * ( hp1 + v * ( hp15 + v * ( hp2 + v * ( hp25 + v * hp3
                ) ) ) ) ) );
        (*hcross)->data->data[idx] = ampfac * v2 * ( hc0 + v * ( hc05
                + v * ( hc1 + v * ( hc15 + v * ( hc2 + v * ( hc25 + v * hc3
                ) ) ) ) ) );

    } /* end loop over time series samples idx */
    return XLAL_SUCCESS;
}
/** @} */

/**
 * Given time series for a binary's orbital dynamical variables,
 * construct the waveform polarizations h+ and hx as a sum of
 * -2 spin-weighted spherical harmonic modes, h_lm.
 * NB: Valid only for non-precessing systems!
 *
 * Implements Equation (11) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
int XLALSimInspiralPNPolarizationWaveformsFromModes(
        REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
               REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
               REAL8TimeSeries *v,       /**< post-Newtonian parameter */
               REAL8TimeSeries *phi,     /**< orbital phase */
               REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
               REAL8 m1,                 /**< mass of companion 1 */
               REAL8 m2,                 /**< mass of companion 2 */
               REAL8 r,                  /**< distance of source */
               REAL8 i,                  /**< inclination of source (rad) */
               int O                     /**< twice post-Newtonian order */
        )
{
    int l, m;
    LAL_CHECK_VALID_SERIES(v, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(phi, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(v, phi, XLAL_FAILURE);
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &v->epoch, 0.0, v->deltaT, &lalStrainUnit, v->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &v->epoch, 0.0, v->deltaT, &lalStrainUnit, v->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(*(*hcross)->data->data));
    for ( l = 2; l <= LAL_PN_MODE_L_MAX; ++l ) {
        for ( m = 1; m <= l; ++m ) {
            COMPLEX16TimeSeries *hmode;
            hmode = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(v, phi, v0, m1, m2, r, O, l, m);
            if ( ! hmode )
                XLAL_ERROR(XLAL_EFUNC);
            if ( XLALSimAddMode(*hplus, *hcross, hmode, i, 0.0, l, m, 1) < 0 )
                XLAL_ERROR(XLAL_EFUNC);
            XLALDestroyCOMPLEX16TimeSeries(hmode);
        }
    }
    return 0;
}

/**
 * Compute the polarizations from all the -2 spin-weighted spherical harmonic
 * modes stored in 'hlms'. Be sure that 'hlms' is the head of the linked list!
 *
 * The computation done is:
 * \f$hp(t) - i hc(t) = \sum_l \sum_m h_lm(t) -2Y_lm(iota,psi)\f$
 *
 * iota and psi are the inclination and polarization angle of the observer
 * relative to the source of GWs.
 */
int XLALSimInspiralPolarizationsFromSphHarmTimeSeries(
    REAL8TimeSeries **hp, /**< Plus polarization time series [returned] */
    REAL8TimeSeries **hc, /**< Cross polarization time series [returned] */
    SphHarmTimeSeries *hlms, /**< Head of linked list of waveform modes */
    REAL8 iota, /**< inclination of viewer to source frame (rad) */
    REAL8 phiRef /**< reference phase (rad) */
    )
{
    int ret;
    SphHarmTimeSeries *ts = hlms;
    size_t length = ts->mode->data->length;
    // Destroy hp, hc TimeSeries if they already exist
    if( (*hp) ) XLALDestroyREAL8TimeSeries( *hp );
    if( (*hc) ) XLALDestroyREAL8TimeSeries( *hc );
    *hp = XLALCreateREAL8TimeSeries("hplus", &(ts->mode->epoch), ts->mode->f0,
                ts->mode->deltaT, &lalStrainUnit, length);
    *hc = XLALCreateREAL8TimeSeries("hplus", &(ts->mode->epoch), ts->mode->f0,
                ts->mode->deltaT, &lalStrainUnit, length);
    memset( (*hp)->data->data, 0, (*hp)->data->length*sizeof(REAL8) );
    memset( (*hc)->data->data, 0, (*hc)->data->length*sizeof(REAL8) );
    while (ts) { // Add the contribution from the current mode to hp, hx...
        // This function adds hlm(t) * Y_lm(incl,phiRef) to (h+ - i hx)(t)
        ret = XLALSimAddMode(*hp, *hc, ts->mode, iota, phiRef, ts->l, ts->m, 0);
        if( ret != XLAL_SUCCESS ) XLAL_ERROR(XLAL_EFUNC);
        ts = ts->next;
    }

    return XLAL_SUCCESS;
}


/**
  Function returning the Fourier domain polarizations for positive frequencies built from the individual modes computed with ChooseFDModes.
	The output should be equivalent to that from ChooseFDWaveform, close to machine precision.
	Some of the AS models use the argument phiRef internally to build the hlms and build the polarizations with an azimuthal angle different from Pi/2 - phiRef.
	In this function we take into account those differences among models.
	For the precessing model IMRPhenomXPHM, since the modes are returned in the J-frame, the polarizations are built using theta_JN instead of the inclination and azimuthal angle = 0.
*/
int XLALSimInspiralPolarizationsFromChooseFDModes(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    const REAL8 m1,                         /**< mass of companion 1 (kg) */
    const REAL8 m2,                         /**< mass of companion 2 (kg) */
    const REAL8 S1x,                        /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                        /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                        /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                        /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                        /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                        /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                   /**< distance of source (m) */
    const REAL8 inclination,                /**< inclination of source (rad) */
    const REAL8 phiRef,                     /**< reference orbital phase (rad) */
    const REAL8 UNUSED longAscNodes,        /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 UNUSED eccentricity,        /**< eccentricity at reference epoch */
    const REAL8 UNUSED meanPerAno,          /**< mean anomaly of periastron */
    // frequency sampling parameters, no default value
    const REAL8 deltaF,                     /**< sampling interval (Hz) */
    const REAL8 f_min,                      /**< starting GW frequency (Hz) */
    const REAL8 f_max,                      /**< ending GW frequency (Hz) */
    REAL8 f_ref,                            /**< Reference frequency (Hz) */
    LALDict *LALparams,                     /**< LAL dictionary containing accessory parameters */
    const Approximant approximant           /**< post-Newtonian approximant to use for waveform production */
    )
{

	/* Adjust the reference frequency for certain precessing approximants:
	 * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
	 * otherwise do nothing */
	FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    int ret = XLAL_SUCCESS;
    REAL8 phiRef_modes = 0;
    REAL8 theta = inclination;
    REAL8 azimuthal = LAL_PI_2 - phiRef;
    REAL8 zeta_polarization = 0;

    switch (approximant)
    {
        case IMRPhenomXHM:
        phiRef_modes = phiRef;
        azimuthal = LAL_PI_2;
        break;

        case IMRPhenomXPHM:
        phiRef_modes = phiRef;
        REAL8 d1=0, d2=0, d3=0, d4=0, d5=0;
        ret = XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame(&d1, &d2, &d3, &theta, &d4, &d5, &zeta_polarization, m1, m2, f_ref, phiRef, inclination, S1x,S1y,S1z, S2x,S2y,S2z, LALparams);
        XLAL_CHECK(XLAL_SUCCESS == ret, XLAL_EFUNC, "Error: XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame failed.\n");
        azimuthal = 0.;

        break;
        case IMRPhenomXO4a:
        phiRef_modes = phiRef;
        d1=0, d2=0, d3=0, d4=0, d5=0;
        ret = XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame(&d1, &d2, &d3, &theta, &d4, &d5, &zeta_polarization, m1, m2, f_ref, phiRef, inclination, S1x,S1y,S1z, S2x,S2y,S2z, LALparams);
        XLAL_CHECK(XLAL_SUCCESS == ret, XLAL_EFUNC, "Error: XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame failed.\n");
        azimuthal = 0.;
        break;
        case SEOBNRv4HM_ROM:

        break;
        case SEOBNRv5_ROM:

        break;
        case SEOBNRv5HM_ROM:

        break;
        case IMRPhenomHM:
        phiRef_modes = phiRef;
        break;
        default:
                XLALPrintError("Approximant not implemented\n");
                XLAL_ERROR(XLAL_EINVAL);
    }


    SphHarmFrequencySeries **hlms = XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms = XLALSimInspiralChooseFDModes(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, deltaF, f_min, f_max, f_ref, phiRef_modes, distance, inclination, LALparams, approximant);
    if(!(*hlms)){
        XLAL_ERROR(XLAL_EFUNC, "Error: XLALSimInspiralChooseFDModes failed\n");
    }

    UINT4 len = (*hlms)->mode->data->length;
    UINT4 offset = 0;

    /* This is to account that ChooseFDModes return the modes for both negative and positve frequencies,
       but here we return the polarizations just for positive frequencies. */
    len = (UINT4) ceil(len/2.);
    offset = len-1;

    /* Create polarizations objects */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("FD hplus",
                    &((*hlms)->mode->epoch), (*hlms)->mode->f0, (*hlms)->mode->deltaF,
                    &((*hlms)->mode->sampleUnits), len);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
    memset((*hptilde)->data->data, 0, (len) * sizeof(COMPLEX16));
    memset((*hctilde)->data->data, 0, (len) * sizeof(COMPLEX16));


    SphHarmFrequencySeries *hlms_tmp = *hlms;

    /* Build the polarizations by summing the modes*/
    while ( hlms_tmp ){
        COMPLEX16 Ylm = XLALSpinWeightedSphericalHarmonic(theta, azimuthal, -2, hlms_tmp->l, hlms_tmp->m);
        COMPLEX16 Ylmstar = conj(Ylm);

        for(UINT4 idx = 0; idx<len; idx++){
            COMPLEX16 hlm = hlms_tmp->mode->data->data[idx + offset];
            COMPLEX16 hlm2 = conj(hlms_tmp->mode->data->data[len - 1 - idx]);
            (*hptilde)->data->data[idx] += 0.5 * (hlm * Ylm + hlm2 * Ylmstar);
            (*hctilde)->data->data[idx] += 0.5 * I * (hlm * Ylm - hlm2 * Ylmstar);
        }
        hlms_tmp = hlms_tmp->next;
    }

    /* Free memory */
    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);
    XLALDestroySphHarmFrequencySeries(hlms_tmp);


    /* Add the correct polarization angle for IMRPhenomXPHM */
    if(fabs(zeta_polarization) > 0)
    {
        COMPLEX16 PhPpolp, PhPpolc;
        REAL8 cosPolFac, sinPolFac;

        cosPolFac = cos(2.0 * zeta_polarization);
        sinPolFac = sin(2.0 * zeta_polarization);

        for (UINT4 i = 0; i < (*hptilde)->data->length; i++)
        {
            PhPpolp = (*hptilde)->data->data[i];
            PhPpolc = (*hctilde)->data->data[i];

            (*hptilde)->data->data[i] = cosPolFac * PhPpolp + sinPolFac * PhPpolc;
            (*hctilde)->data->data[i] = cosPolFac * PhPpolc - sinPolFac * PhPpolp;
        }
    }


  /* This final rotation is taken from ChooseFDWaveform */
    REAL8 polariz=longAscNodes;
    if (polariz) {
        COMPLEX16 tmpP,tmpC;
        for (UINT4 idx=0;idx<(*hptilde)->data->length;idx++) {
tmpP=(*hptilde)->data->data[idx];
tmpC=(*hctilde)->data->data[idx];
(*hptilde)->data->data[idx] =cos(2.*polariz)*tmpP+sin(2.*polariz)*tmpC;
(*hctilde)->data->data[idx]=cos(2.*polariz)*tmpC-sin(2.*polariz)*tmpP;
        }
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
    if (XLALSimInspiralWaveformParamsLookupEnableLIV(LALparams))
        ret = XLALSimLorentzInvarianceViolationTerm(hptilde, hctilde, m1/LAL_MSUN_SI, m2/LAL_MSUN_SI, distance, LALparams);
    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

    return ret;
}



/**
    Return polarizations for positive frequencies built by summing the individual modes present
    in the input array SphHarmFrequencySeries *hlms	computed with ChooseFDModes.
    Notice that in general the output may not be close to machine precision with ChooseFDWaveform due to differences in
    computing the h_lms and in the use of the azimuthal angle. 
    For IMRPhenomXPHM, the argument theta should not be the inclination but theta_JN and phiRef should be 0.
*/
int XLALSimInspiralPolarizationsFromSphHarmFrequencySeries(
    COMPLEX16FrequencySeries **hp, /**< Plus polarization frequency series [returned] */
    COMPLEX16FrequencySeries **hc, /**< Cross polarization frequency series [returned] */
    SphHarmFrequencySeries *hlms,  /**< Head of linked list of waveform modes */
    REAL8 theta, 									 /**< Polar angle for the Ylms (rad): Inclination for h_lms in L0-frame and theta_JN for J-frame. */
    REAL8 phi 										 /**< Azimuthal angle for the Ylms (rad): pi/2 - phiRef for AS models and 0 for precessing. */
    )
{
    if (!(hlms)){ XLAL_ERROR(XLAL_EFUNC, "SphHarmFrequencySeires object empty.\n");}

    SphHarmFrequencySeries *fs = hlms;
    UINT4 len = fs->mode->data->length;
    INT4 offset = 0;
    /* This is to account that ChooseFDModes return the modes for both negative and positve frequencies,
        but here we return the polarizations just for positive frequencies. */
    len = (UINT4) ceil(len/2.);
    offset = len-1;
    // Destroy hp, hc FrequencySeries if they already exist
    if( (*hp) ) XLALDestroyCOMPLEX16FrequencySeries( *hp );
    if( (*hc) ) XLALDestroyCOMPLEX16FrequencySeries( *hc );
    *hp = XLALCreateCOMPLEX16FrequencySeries("hplus", &(fs->mode->epoch), fs->mode->f0,
                fs->mode->deltaF, &(fs->mode->sampleUnits), len);
    *hc = XLALCreateCOMPLEX16FrequencySeries("hcross", &(fs->mode->epoch), fs->mode->f0,
                fs->mode->deltaF, &(fs->mode->sampleUnits), len);
    memset( (*hp)->data->data, 0, (*hp)->data->length*sizeof(COMPLEX16) );
    memset( (*hc)->data->data, 0, (*hc)->data->length*sizeof(COMPLEX16) );


    /* Build the polarizations by summing the modes */
    while ( fs ){
        COMPLEX16 Ylm = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, fs->l, fs->m);
        COMPLEX16 Ylmstar = conj(Ylm);

        for(UINT4 idx = 0; idx<len; idx++){
            COMPLEX16 hlm = fs->mode->data->data[idx + offset];
            COMPLEX16 hlm2 = conj(fs->mode->data->data[len - 1 - idx]);
            (*hp)->data->data[idx] += 0.5 * (hlm * Ylm + hlm2 * Ylmstar);
            (*hc)->data->data[idx] += 0.5 * I * (hlm * Ylm - hlm2 * Ylmstar);
        }
        fs = fs->next;
    }
    XLALDestroySphHarmFrequencySeries(fs);

    return XLAL_SUCCESS;
}


/**
 * Given time series for a binary's orbital dynamical variables,
 * computes the radial and angular orbital motion; then constructs
 * the waveform polarizations h+ and hx in terms of the relative
 * separation r, the `true anomaly` phi and their time derivatives.
 * NB: Valid only for non-spinning binaries in inspiralling!
 * and precessing eccentric orbits!
 *
 * Implements Equations (3.7a) - (3.7c) and Equations (B2a) - (B2d), (B4a), (B4b) of:
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * https://dcc.ligo.org/P1500148-v1
 *
 * (note that above equations use x = v^2 as the PN expansion parameter)
 *
 * as well as Equations (6a) and (6b) of:
 * Thibault Damour, Achamveedu Gopakumar, and Bala R. Iyer,
 * \"Phasing of gravitational waves from inspiralling eccentric
 * binaries\", Phys. Rev. D 70 064028 (2004);
 * arXiv:gr-qc/0404128.
 */
int XLALSimInspiralPNPolarizationWaveformsEccentric(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
    REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
    REAL8TimeSeries *V,       /**< post-Newtonian (PN) parameter */
    REAL8TimeSeries *Ecc,     /**< orbital eccentricity */
    REAL8TimeSeries *U,       /**< eccentric anomaly of quasi-Keplerian orbit */
    REAL8TimeSeries *Phi,     /**< orbital phase */
    REAL8 m1,                 /**< mass of companion 1 (kg) */
    REAL8 m2,                 /**< mass of companion 2 (kg) */
    REAL8 r,                  /**< distance of source (m) */
    REAL8 i,                  /**< inclination of source (rad) */
    int ampO,                 /**< twice PN order of the amplitude */
    int ph_O                  /**< twice PN order of the phase */
    )
{
    REAL8 mt, eta, dist, ampfac, phi, v, et, u;
    REAL8 dt, OTS;
    REAL8 CF_rdot, CF_rphidot, CF_Z, r_sc, rdot_sc, phidot;
    REAL8 rdot, rphidot, Z;
    REAL8 hp0;
    REAL8 hc0;
    REAL8 ci, si, ci2, si2;
    INT4 idx, len;

    /* Sanity check input time series */
    LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Ecc, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(U, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Ecc, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, U, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length
            * sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length
            * sizeof(*(*hcross)->data->data));

    mt = (m1 + m2)/LAL_MSUN_SI;
    eta = m1 * m2 / pow(m1 + m2, 2.0); // symmetric mass ratio - '\nu' in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = mt * LAL_MTSUN_SI * eta / dist;

    /*
     * cosines and sines of inclination between
     * line of sight (N) and binary orbital angular momentum (L_N)
     */
    ci = cos(i);  	si = sin(i);
    ci2 = ci*ci;
    si2 = si*si;

    /* loop over time steps and compute first radial and angular orbital variables,
     * then the waveform polarizations h+ and hx */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi = Phi->data->data[idx]; 	et = Ecc->data->data[idx];	u = U->data->data[idx];	v = V->data->data[idx];

    /* Some abbreviated functions in terms of u and et. */
    dt	= 1. - et * cos(u);
    OTS	= sqrt(1. - et * et);

    /*
         * First set the functions for the dimensionless orbital variables
     * (1/c)*dr/dt, (r/c)*dphi/dt, Z = G*m/(r*c^2) to 0.
     * Then use a switch to set proper non-zero values
         * up to order ph_O for the contributions to the dimensionless variables. Note we
         * fall through the PN orders and break only after Newt. order
         */

    rdot = rphidot = Z = 0;

    /* Common factors in PN accurate expressions for orbital variables */
    CF_rdot = et * v*sin(u)/dt;
    CF_rphidot = OTS * v / dt;
    CF_Z = pow(v, 2.0) / dt;

    switch( ph_O )
        {
        /* case LAL_PNORDER_FOUR: */
            case 8:
                XLALPrintError("XLAL Error - %s: dynamical variables not known "
                        "to PN order %d\n", __func__, ph_O );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_THREE_POINT_FIVE: */
            case 7:
                XLALPrintError("XLAL Error - %s: dynamical variables not known "
                        "to PN order %d\n", __func__, ph_O );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_THREE: */
            case 6:
        XLALPrintError("XLAL Error - %s: dynamical variables not yet implemented "
                        "to PN order %d\n", __func__, ph_O );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_TWO_POINT_FIVE: */
            case 5:
        XLALPrintError("XLAL Error - %s: dynamical variables not yet implemented "
                        "to PN order %d\n", __func__, ph_O );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        case -1: // Highest known PN order - move if higher terms added!
            /* case LAL_PNORDER_TWO: */
            case 4:
        r_sc 	= 1.
              //==============================================================================
              + ((-24. + dt * (18. - 7. * eta) + 9. * eta + pow(et, 2.0) * (24. - 9. * eta
              + dt * (-6. + 7. * eta))) * pow(v, 2.0)) / (6. * dt * pow(OTS, 2.0))
              //==============================================================================
              + ((-288. + 765. * eta - 27. * pow(eta,  2.0)
              + pow(et, 4.0) * (261. * eta - 27 * pow(eta, 2.0))
              + pow(et, 2.0) * (288. - 1026. * eta + 54. * pow(eta, 2.0))
              + (-540. + pow(et, 2.0) * (540. - 216. * eta)
              + 216. * eta) * OTS + dt * (648. - 567. * eta + 35. * pow(eta, 2.0)
              + pow(et, 2.0) * (468. + 150. * eta - 70. * pow(eta, 2.0))
              + pow(et, 4.0) * (72. - 231. * eta + 35. * pow(eta, 2.0))
              + (180. - 72. * eta + pow(et, 2.0) * (-180. + 72. * eta)) * OTS)) * pow(v, 4.0))
              / (72. * dt * pow(OTS, 4.0));
              //==============================================================================
        phidot = 1.
              //===========================================================================================
              + (-1. + dt + pow(et, 2.0)) * (-4. + eta) * pow(v, 2.0) / (dt * pow(OTS, 2.0))
              //===========================================================================================
              + (-6. * pow(1. - pow(et, 2.0), 3.0) * eta * (3. + 2. * eta)
              + pow(dt, 3.0) * (42. + 22. * eta + 8.* pow(eta, 2.0)
              + pow(et, 2.0) * (-147. + 8. * eta - 14. * pow(eta, 2.0)))
              + dt * (108. + 63. * eta + 33. * pow(eta, 2.0)
              + pow(et, 2.0) * (-216. - 126. * eta - 66. * pow(eta, 2.0))
              + pow(et, 4.0) * (108. + 63. * eta + 33. * pow(eta, 2.0)))
              + pow(dt, 2.0) * (-240. - 31. * eta - 29. * pow(eta, 2.0)
              + pow(et, 4.0) * (-48. + 17. * eta - 17. * pow(eta, 2.0))
              + pow(et, 2.0) * (288. + 14. * eta + 46. * pow(eta,2)))
              + 18. * pow(dt, 2.0) * (-2. + dt + 2. * pow(et, 2.0)) * (-5. + 2. * eta) * OTS) * pow(v, 4.0)
              / (12. * pow(dt, 3.0) * pow(OTS, 4.0));
              //===========================================================================================
        rdot_sc = 1.
              //==========================================================================================
              + (-7. * eta + pow(et, 2.0) * (-6. + 7. * eta)) * pow(v, 2.0) / (6. * pow(OTS, 2.0))
              //==========================================================================================
              + (-135. * eta + 9. * pow(eta, 2.0) + pow(et, 2.0) * (405. * eta - 27. * pow(eta, 2.0))
              + pow(et, 6.0) * (135. * eta - 9. * pow(eta, 2.0))
              + pow(et, 4.0) * (-405. * eta + 27 * pow(eta, 2.0))
              + dt * (-540. + 351. * eta - 9. * pow(eta, 2.0)
              + pow(et, 4.0) * (-540. + 351. * eta - 9. * pow(eta, 2.0))
              + pow(et, 2.0) * (1080. - 702. * eta + 18. * pow(eta, 2.0)))
              + pow(dt, 3.0) * (-324. + 189. * eta + 35. * pow(eta, 2.0)
              + pow(et, 2.0) * (-234. + 366. * eta - 70. * pow(eta, 2.0))
              + pow(et, 4.0) * (72. - 231. * eta + 35. * pow(eta, 2.0)))
              - 36 * pow(dt, 2.0) * (3. + dt) * (1 - pow(et, 2.0)) * (-5. + 2. * eta) * OTS) * pow(v, 4.0)
              / (72. * pow(dt, 3.0) * pow(OTS, 4.0));
              //==========================================================================================
                break;
            /* case LAL_PNORDER_ONE_POINT_FIVE: */
            case 3:
        r_sc 	= 1.
              //===========================================================================
              + ((-24. + dt * (18. - 7. * eta) + 9. * eta + pow(et, 2.0) * (24. - 9. * eta
              + dt * (-6. + 7. * eta))) * pow(v, 2.0)) / (6. * dt * pow(OTS, 2.0));
              //===========================================================================
        phidot	= 1.
              //==============================================================================
              + (-1. + dt + pow(et, 2.0)) * (-4. + eta) * pow(v, 2.0) / (dt * pow(OTS, 2.0));
              //==============================================================================
        rdot_sc = 1.
              //====================================================================================
              + (-7. * eta + pow(et, 2.0) * (-6. + 7. * eta)) * pow(v, 2.0) / (6. * pow(OTS, 2.0));
              //====================================================================================
                break;
            /* case LAL_PNORDER_ONE: */
            case 2:
        r_sc 	= 1.
              //===========================================================================
              + ((-24. + dt * (18. - 7. * eta) + 9. * eta + pow(et, 2.0) * (24. - 9. * eta
              + dt * (-6. + 7. * eta))) * pow(v, 2.0)) / (6. * dt * pow(OTS, 2.0));
              //===========================================================================
        phidot	= 1.
              //==============================================================================
              + (-1. + dt + pow(et, 2.0)) * (-4. + eta) * pow(v, 2.0) / (dt * pow(OTS, 2.0));
              //==============================================================================
        rdot_sc = 1.
              //====================================================================================
              + (-7. * eta + pow(et, 2.0) * (-6. + 7. * eta)) * pow(v, 2.0) / (6. * pow(OTS, 2.0));
              //====================================================================================
                break;
            /*case LAL_PNORDER_HALF:*/
            case 1:
        r_sc	= 1.;
        phidot	= 1.;
        rdot_sc = 1.;
                break;
        /*case LAL_PNORDER_NEWTONIAN:*/
            case 0:
                r_sc	= 1.;
        phidot	= 1.;
        rdot_sc = 1.;
                break;
            default:
                XLALPrintError("XLAL Error - %s: Invalid phase PN order %d\n",
                        __func__, ph_O );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        }

        /* Dimensionless rdot/c, (r/c)*phidot, Z = G*m/(r*c^2) entering the h+/x coefficients*/
        rdot = CF_rdot * rdot_sc;
    rphidot = CF_rphidot * r_sc * phidot;
        Z = CF_Z/r_sc;

        /*
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hp0 = 0.;
        hc0 = 0.;

        switch( ampO )
        {
            /* case LAL_PNORDER_THREE_POINT_FIVE: */
            case 7:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_THREE: */
            case 6:
        XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_TWO_POINT_FIVE: */
            case 5:
        XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_TWO: */
            case 4:
        XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_ONE_POINT_FIVE: */
            case 3:
        XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /* case LAL_PNORDER_ONE: */
            case 2:
        XLALPrintError("XLAL Error - %s: Amp. corrections not yet implemented "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            /*case LAL_PNORDER_HALF:*/
            case 1:
        XLALPrintError("XLAL Error - %s: Amp. corrections not yet implemented "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        case -1: // Highest known PN order - move if higher terms added!
        /*case LAL_PNORDER_NEWTONIAN:*/
            case 0:
          hp0 = - (si2 * (-pow(rphidot, 2.0) - pow(rdot, 2.0) + Z))
            - (1. + ci2) * ((pow(rphidot, 2.0) - pow(rdot, 2.0) + Z) * cos(2. * phi)
            + (2. * rphidot * rdot) * sin(2. * phi));
          hc0 = - 2. * ci * (-2. * rphidot * rdot * cos(2. * phi)
              + (pow(rphidot, 2.0) - pow(rdot, 2.0) + Z) * sin(2. * phi));
                break;
            default:
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d\n",
                        __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
      (*hplus)->data->data[idx] = ampfac*hp0;
      (*hcross)->data->data[idx] = ampfac*hc0;

    } /* end loop over time series samples idx */
    return XLAL_SUCCESS;
}

/**
 * Computes polarizations h+ and hx for a spinning, precessing binary
 * when provided time series of all the dynamical quantities.
 * Amplitude can be chosen between 1.5PN and Newtonian orders (inclusive).
 *
 * Based on K.G. Arun, Alesssandra Buonanno, Guillaume Faye and Evan Ochsner
 * \"Higher-order spin effects in the amplitude and phase of gravitational
 * waveforms emitted by inspiraling compact binaries: Ready-to-use
 * gravitational waveforms\", Phys Rev. D 79, 104023 (2009), arXiv:0810.5336
 *
 * HOWEVER, the formulae have been adapted to use the output of the so-called
 * \"Frameless\" convention for evolving precessing binary dynamics,
 * which is not susceptible to hitting coordinate singularities.
 *
 * FIXME: Clean up and commit Mathematica NB Showing correctness. Cite here.
 *
 * NOTE: The vectors MUST be given in the so-called radiation frame where
 * Z is the direction of propagation, X is the principal '+' axis and Y = Z x X
 * For different convention (Z is the direction of initial total angular
 * momentum, useful for GRB and comparison to NR, see XLALSimSpinInspiralGenerator())
 */
int XLALSimInspiralPrecessingPolarizationWaveforms(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
    REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
    REAL8TimeSeries *V,       /**< post-Newtonian parameter */
    REAL8TimeSeries *Phi,     /**< orbital phase */
    REAL8TimeSeries *S1x,	  /**< Spin1 vector x component */
    REAL8TimeSeries *S1y,	  /**< Spin1 vector y component */
    REAL8TimeSeries *S1z,	  /**< Spin1 vector z component */
    REAL8TimeSeries *S2x,	  /**< Spin2 vector x component */
    REAL8TimeSeries *S2y,	  /**< Spin2 vector y component */
    REAL8TimeSeries *S2z,	  /**< Spin2 vector z component */
    REAL8TimeSeries *LNhatx,  /**< unit orbital ang. mom. x comp. */
    REAL8TimeSeries *LNhaty,  /**< unit orbital ang. mom. y comp. */
    REAL8TimeSeries *LNhatz,  /**< unit orbital ang. mom. z comp. */
    REAL8TimeSeries *E1x,	  /**< orbital plane basis vector x comp. */
    REAL8TimeSeries *E1y,	  /**< orbital plane basis vector y comp. */
    REAL8TimeSeries *E1z,	  /**< orbital plane basis vector z comp. */
    REAL8 m1,                 /**< mass of companion 1 (kg) */
    REAL8 m2,                 /**< mass of companion 2 (kg) */
    REAL8 r,                  /**< distance of source (m) */
    INT4 ampO	 	  /**< twice amp. post-Newtonian order */
    )
{
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z, lnhx, lnhy, lnhz;
    REAL8 e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, lx, ly, lz;
    REAL8 nx2, ny2, nz2, nz3, lx2, ly2, lz2, lz3;
    REAL8 hplus0, hcross0, hplus05, hcross05, hplus1, hcross1;
    REAL8 hplus15, hcross15, hplusSpin1, hcrossSpin1;
    REAL8 hplusSpin15, hcrossSpin15, hplusTail15, hcrossTail15;
    REAL8 M, eta, dm, phi, v, v2, dist, ampfac;
    INT4 idx, len;

    /* Macros to check time series vectors */
    if ( !(hplus) || !(hcross) ) {
      XLALPrintError("** XLALSimInspiralPrecessingPolarizationWaveform error: **h+ and **hx are expected to be non NULL!\n");
      XLAL_ERROR( XLAL_FAILURE );
    }
    LAL_CHECK_VALID_SERIES(V, 			XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1z, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2z, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatx, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhaty, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatz, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1z, 		XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatx, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhaty, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatz, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1z, 	XLAL_FAILURE);

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0,
            (*hplus)->data->length*sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0,
            (*hcross)->data->length*sizeof(*(*hcross)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    dm = (m1 - m2) / M;    // frac. mass difference - \delta m/m in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;

    /* loop over time steps and compute polarizations h+ and hx */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi  = Phi->data->data[idx]; 	v = V->data->data[idx];     v2 = v * v;
        lnhx = LNhatx->data->data[idx]; e1x = E1x->data->data[idx];
        lnhy = LNhaty->data->data[idx];	e1y = E1y->data->data[idx];
        lnhz = LNhatz->data->data[idx];	e1z = E1z->data->data[idx];
        s1x  = S1x->data->data[idx];	s2x = S2x->data->data[idx];
        s1y  = S1y->data->data[idx];	s2y = S2y->data->data[idx];
        s1z  = S1z->data->data[idx];	s2z = S2z->data->data[idx];

        /* E2 = LNhat x E1 */
        e2x = lnhy*e1z - lnhz*e1y;
        e2y = lnhz*e1x - lnhx*e1z;
        e2z = lnhx*e1y - lnhy*e1x;

        /* Unit orbital separation vector */
        nx = e1x*cos(phi) + e2x*sin(phi);
        ny = e1y*cos(phi) + e2y*sin(phi);
        nz = e1z*cos(phi) + e2z*sin(phi);

        /* Unit inst. orbital velocity vector */
        lx = e2x*cos(phi) - e1x*sin(phi);
        ly = e2y*cos(phi) - e1y*sin(phi);
        lz = e2z*cos(phi) - e1z*sin(phi);

        /* Powers of vector components */
        nx2 = nx*nx;	ny2 = ny*ny;	nz2 = nz*nz;	nz3 = nz*nz2;
        lx2 = lx*lx;	ly2 = ly*ly;	lz2 = lz*lz;	lz3 = lz*lz2;

        /*
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hplus0 = hplus05 = hplus1 = hplus15 = hplusTail15 = 0.;
        hcross0 = hcross05 = hcross1 = hcross15 = hcrossTail15 = 0.;
        hplusSpin1 = hplusSpin15 = hcrossSpin1 = hcrossSpin15 = 0.;

        switch( ampO )
        {
            /*
             * case LAL_PNORDER_THREE_POINT_FIVE:
             * case LAL_PNORDER_THREE:
             * case LAL_PNORDER_TWO_POINT_FIVE:
             * case LAL_PNORDER_TWO:
             */
            case 7:
            case 6:
            case 5:
            case 4:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d, highest is %d\n", __func__, ampO,
                        MAX_PRECESSING_AMP_PN_ORDER );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            case -1: /* Use highest known PN order - move if new orders added */
            /*case LAL_PNORDER_ONE_POINT_FIVE:*/
            case 3:
                /* 1.5PN non-spinning amp. corrections */
                hplus15 = (dm*(2*lx*nx*nz*(-95 + 90*lz2 - 65*nz2
                        - 2*eta*(-9 + 90*lz2 - 65*nz2)) - 2*ly*ny*nz
                        * (-95 + 90*lz2 - 65*nz2 - 2*eta*(-9 + 90*lz2 - 65*nz2))
                        + 6*lx2*lz*(13 - 4*lz2 + 29*nz2 + eta*(-2 + 8*lz2
                        - 58*nz2)) - 6*ly2*lz*(13 - 4*lz2 + 29*nz2 + eta
                        * (-2 + 8*lz2 - 58*nz2)) - lz*(nx2 - ny2)*(83 - 6*lz2
                        + 111*nz2 + 6*eta*(-1 + 2*lz2 - 37*nz2))))/24.;
                hcross15 = (dm*(lz*(6*(19 - 4*eta)*lx*ly + (-101 + 12*eta)
                        * nx*ny) + (-149 + 36*eta) * (ly*nx + lx*ny)*nz
                        + 6*(-3 + eta) * (2*lx*ly*lz - lz*nx*ny - 3*ly*nx*nz
                        - 3*lx*ny*nz) + (1 - 2*eta) * (6*lz3*(-4*lx*ly + nx*ny)
                        + 90*lz2*(ly*nx + lx*ny)*nz + 3*lz*(58*lx*ly
                        - 37*nx*ny)*nz2 - 65*(ly*nx + lx*ny)*nz3)))/12.;
                /* 1.5PN spinning amp. corrections */
                hplusSpin15 = (6*lz*ny*s1x + 6*dm*lz*ny*s1x - 3*eta*lz*ny*s1x
                        + 2*ly2*lnhy*s1y + 2*dm*ly2*lnhy*s1y
                        + 2*eta*ly2*lnhy*s1y + 6*lz*nx*s1y + 6*dm*lz*nx*s1y
                        - 3*eta*lz*nx*s1y + 8*lnhy*nx2*s1y + 8*dm*lnhy*nx2*s1y
                        - eta*lnhy*nx2*s1y - 8*lnhy*ny2*s1y - 8*dm*lnhy*ny2*s1y
                        + eta*lnhy*ny2*s1y + 2*ly2*lnhz*s1z + 2*dm*ly2*lnhz*s1z
                        + 2*eta*ly2*lnhz*s1z - 6*ly*nx*s1z - 6*dm*ly*nx*s1z
                        - 9*eta*ly*nx*s1z + 8*lnhz*nx2*s1z + 8*dm*lnhz*nx2*s1z
                        - eta*lnhz*nx2*s1z - 8*lnhz*ny2*s1z - 8*dm*lnhz*ny2*s1z
                        + eta*lnhz*ny2*s1z + 6*lz*ny*s2x - 6*dm*lz*ny*s2x
                        - 3*eta*lz*ny*s2x + lnhx*(2*ly2*((1 + dm + eta)*s1x
                        + (1 - dm + eta)*s2x) + (nx2 - ny2)*((8 + 8*dm - eta)
                        * s1x - (-8 + 8*dm + eta)*s2x)) + 2*ly2*lnhy*s2y
                        - 2*dm*ly2*lnhy*s2y + 2*eta*ly2*lnhy*s2y + 6*lz*nx*s2y
                        - 6*dm*lz*nx*s2y - 3*eta*lz*nx*s2y + 8*lnhy*nx2*s2y
                        - 8*dm*lnhy*nx2*s2y - eta*lnhy*nx2*s2y - 8*lnhy*ny2*s2y
                        + 8*dm*lnhy*ny2*s2y + eta*lnhy*ny2*s2y + 2*ly2*lnhz*s2z
                        - 2*dm*ly2*lnhz*s2z + 2*eta*ly2*lnhz*s2z - 6*ly*nx*s2z
                        + 6*dm*ly*nx*s2z - 9*eta*ly*nx*s2z + 8*lnhz*nx2*s2z
                        - 8*dm*lnhz*nx2*s2z - eta*lnhz*nx2*s2z - 8*lnhz*ny2*s2z
                        + 8*dm*lnhz*ny2*s2z + eta*lnhz*ny2*s2z - 3*lx*ny
                        * ((2 + 2*dm + 3*eta)*s1z + (2 - 2*dm + 3*eta)*s2z)
                        - 2*lx2*(lnhx*((1 + dm + eta)*s1x + (1 - dm + eta)*s2x)
                        + lnhy*((1 + dm + eta)*s1y + (1 - dm + eta)*s2y)
                        + lnhz*((1 + dm + eta)*s1z + (1 - dm + eta)*s2z)))/3.;
                hcrossSpin15 = (-3*lz*(nx*((2 + 2*dm - eta)*s1x
                        - (-2 + 2*dm + eta)*s2x) + ny*((-2 - 2*dm + eta)*s1y
                        + (-2 + 2*dm + eta)*s2y)) + ny*(-6*ly*s1z - 6*dm*ly*s1z
                        - 9*eta*ly*s1z + 16*lnhz*nx*s1z + 16*dm*lnhz*nx*s1z
                        - 2*eta*lnhz*nx*s1z + 2*lnhx*nx*((8 + 8*dm - eta)*s1x
                        - (-8 + 8*dm + eta)*s2x) + 2*lnhy*nx*((8 + 8*dm - eta)
                        * s1y - (-8 + 8*dm + eta)*s2y) - 6*ly*s2z + 6*dm*ly*s2z
                        - 9*eta*ly*s2z + 16*lnhz*nx*s2z - 16*dm*lnhz*nx*s2z
                        - 2*eta*lnhz*nx*s2z) - lx*(4*lnhx*ly*((1 + dm + eta)*s1x
                        + (1 - dm + eta)*s2x) - 3*nx*((2 + 2*dm + 3*eta)*s1z
                        + (2 - 2*dm + 3*eta)*s2z) + 4*ly*(lnhy*((1 + dm + eta)
                        * s1y + (1 - dm + eta)*s2y) + lnhz*((1 + dm + eta)*s1z
                        + (1 - dm + eta)*s2z))))/3.;
                /* 1.5PN tail amp. corrections */
                hplusTail15 = 2*((lx2 - ly2 - nx2 + ny2)*LAL_PI);
                hcrossTail15 = 4*((lx*ly - nx*ny)*LAL_PI);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif

            /*case LAL_PNORDER_ONE:*/
            case 2:
                /* 1PN non-spinning amp. corrections */
                hplus1 = (-13*lx2 + 13*ly2 + 6*lx2*lz2 - 6*ly2*lz2
                        + 13*(nx2 - ny2) - 2*lz2*(nx2 - ny2) - 32*lx*lz*nx*nz
                        + 32*ly*lz*ny*nz - 14*lx2*nz2 + 14*ly2*nz2
                        + 10*(nx2 - ny2)*nz2)/6. + (eta*(lx2 - 18*lx2*lz2
                        + 96*lx*lz*nx*nz - 96*ly*lz*ny*nz + 42*lx2*nz2
                        + ly2*(-1 + 18*lz2 - 42*nz2) + (nx2 - ny2)
                        * (-1 + 6*lz2 - 30*nz2)))/6.;
                hcross1 = (eta*(lx*ly - nx*ny - 6*(lz2*(3*lx*ly - nx*ny)
                        - 8*lz*(ly*nx + lx*ny)*nz + (-7*lx*ly
                        + 5*nx*ny)*nz2)))/3. + (-13*(lx*ly - nx*ny)
                        + 2*(lz2*(3*lx*ly - nx*ny) - 8*lz*(ly*nx + lx*ny)*nz
                        + (-7*lx*ly + 5*nx*ny)*nz2))/3.;
                /* 1PN spinning amp. corrections */
                hplusSpin1 = (-(ny*((1 + dm)*s1x + (-1 + dm)*s2x))
                        - nx*((1 + dm)*s1y + (-1 + dm)*s2y))/2.;
                hcrossSpin1 = (nx*((1 + dm)*s1x + (-1 + dm)*s2x)
                        - ny*((1 + dm)*s1y + (-1 + dm)*s2y))/2.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif

            /*case LAL_PNORDER_HALF:*/
            case 1:
                /* 0.5PN non-spinning amp. corrections */
                hplus05 = (dm*(-2*lx2*lz + 2*ly2*lz + lz*(nx2 - ny2)
                        + 6*lx*nx*nz - 6*ly*ny*nz))/2.;
                hcross05 = dm*(-2*lx*ly*lz + lz*nx*ny
                    + 3*ly*nx*nz + 3*lx*ny*nz);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            /*case LAL_PNORDER_NEWTONIAN:*/
            case 0:
                /* Newtonian order polarizations */
                hplus0 = lx2 - ly2 - nx2 + ny2;
                hcross0 = 2*lx*ly - 2*nx*ny;

                break;
            default:
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d\n",
                        __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
        (*hplus)->data->data[idx] = ampfac * v2 * ( hplus0
                + v * ( hplus05 + v * ( hplus1 + hplusSpin1
                + v * ( hplus15 + hplusSpin15 + hplusTail15 ) ) ) );
        (*hcross)->data->data[idx] = ampfac * v2 * ( hcross0
                + v * ( hcross05 + v * ( hcross1 + hcrossSpin1
                + v * ( hcross15 + hcrossSpin15 + hcrossTail15 ) ) ) );
    } /* end of loop over time samples, idx */
    return XLAL_SUCCESS;
}

/**
 * Computes polarizations h+ and hx for a spinning, precessing binary
 * when provided a single value of all the dynamical quantities.
 * Amplitude can be chosen between 1.5PN and Newtonian orders (inclusive).
 *
 * Based on K.G. Arun, Alesssandra Buonanno, Guillaume Faye and Evan Ochsner
 * \"Higher-order spin effects in the amplitude and phase of gravitational
 * waveforms emitted by inspiraling compact binaries: Ready-to-use
 * gravitational waveforms\", Phys Rev. D 79, 104023 (2009), arXiv:0810.5336
 *
 * HOWEVER, the formulae have been adapted to use the output of the so-called
 * \"Frameless\" convention for evolving precessing binary dynamics,
 * which is not susceptible to hitting coordinate singularities.
 *
 * This has been written to reproduce XLALSimInspiralPrecessingPolarizationWaveforms.
 * If hplus and hcross are the output of XLALSimInspiralPrecessingPolarizationWaveforms,
 * and hp(n) and hc(n) the output of this function for a given harmonic number, then
 *
 * hplus = sum_{n=0}^5 hp(n)*exp(-i*n*Phi) + c.c.
 * hcross = sum_{n=0}^5 hc(n)*exp(-i*n*Phi) + c.c.
 *
 * NOTE: The vectors MUST be given in the so-called radiation frame where
 * Z is the direction of propagation, X is the principal '+' axis and Y = Z x X
 * For different convention (Z is the direction of initial total angular
 * momentum, useful for GRB and comparison to NR, see XLALSimSpinInspiralGenerator())
 * FIXME take out v0 as it can be absorbed in a 4PN additional phase term
 * see discussion in sec. VIII of  Class.Quant.Grav. 25 (2008) 165003,
 * arXiv 0802.1249.
 */
int XLALSimInspiralPrecessingPolarizationWaveformHarmonic(
        COMPLEX16 *hplus,  /**< +-polarization waveform [returned] */
        COMPLEX16 *hcross, /**< x-polarization waveform [returned] */
        REAL8 v,       /**< post-Newtonian parameter */
        REAL8 s1x,     /**< Spin1 vector x component */
        REAL8 s1y,     /**< Spin1 vector y component */
        REAL8 s1z,     /**< Spin1 vector z component */
        REAL8 s2x,     /**< Spin2 vector x component */
        REAL8 s2y,     /**< Spin2 vector y component */
        REAL8 s2z,     /**< Spin2 vector z component */
        REAL8 lnhx,  /**< unit orbital ang. mom. x comp. */
        REAL8 lnhy,  /**< unit orbital ang. mom. y comp. */
        REAL8 lnhz,  /**< unit orbital ang. mom. z comp. */
        REAL8 e1x,     /**< orbital plane basis vector x comp. */
        REAL8 e1y,     /**< orbital plane basis vector y comp. */
        REAL8 e1z,     /**< orbital plane basis vector z comp. */
        REAL8 dm,                 /**< dimensionless mass difference (m1 - m2)/(m1 + m2) > 0 */
        REAL8 eta,                /**< symmetric mass ratio m1*m2/(m1 + m2)^2 */
        REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
        INT4 n,                   /**< harmonic number */
        INT4 ampO                 /**< twice amp. post-Newtonian order */
        )
{
  /* E2 = LNhat x E1 */
  REAL8 e2x = lnhy*e1z - lnhz*e1y;
  REAL8 e2y = lnhz*e1x - lnhx*e1z;
  REAL8 e2z = lnhx*e1y - lnhy*e1x;

  REAL8 v2 = v*v;
  REAL8 v3 = v2*v;
  REAL8 v4 = v3*v;
  REAL8 v5 = v4*v;

  REAL8 twom1 = (1. + dm);
  REAL8 twom2 = (1. - dm);

  REAL8 a1x = s1x*twom1;
  REAL8 a1y = s1y*twom1;
  REAL8 a1z = s1z*twom1;
  REAL8 a2x = s2x*twom2;
  REAL8 a2y = s2y*twom2;
  REAL8 a2z = s2z*twom2;

  REAL8 logfac;

  /*
   * First set all h+/x coefficients to 0. Then use a switch to
   * set proper non-zero values up to order ampO. Note we
   * fall through the PN orders and break only after Newt. order
   */
  *hplus = 0.;
  *hcross = 0.;

  REAL8 fact1, fact2, fact3, fact4, fact5, fact6, fact7, fact8, fact9;

  REAL8 e1xe1x = e1x*e1x;
  REAL8 e1xe1y = e1x*e1y;
  REAL8 e1xe1z = e1x*e1z;
  REAL8 e1ye1y = e1y*e1y;
  REAL8 e1ye1z = e1y*e1z;
  REAL8 e1ze1z = e1z*e1z;

  REAL8 e2xe2x = e2x*e2x;
  REAL8 e2xe2y = e2x*e2y;
  REAL8 e2xe2z = e2x*e2z;
  REAL8 e2ye2y = e2y*e2y;
  REAL8 e2ye2z = e2y*e2z;
  REAL8 e2ze2z = e2z*e2z;

  REAL8 e1xe2x = e1x*e2x;
  REAL8 e1xe2y = e1x*e2y;
  REAL8 e1ye2x = e1y*e2x;
  REAL8 e1xe2z = e1x*e2z;
  REAL8 e1ze2x = e1z*e2x;
  REAL8 e1ye2y = e1y*e2y;
  REAL8 e1ye2z = e1y*e2z;
  REAL8 e1ze2y = e1z*e2y;
  REAL8 e1ze2z = e1z*e2z;

  switch(n) // harmonic number
  {
    case 0:
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
        case 2:
        case 1:
          fact1 = v3*0.125;
          fact2 = 7. + dm;
          fact3 = 7. - dm;
          fact4 = a1x*fact2 + a2x*fact3;
          fact5 = a1y*fact2 + a2y*fact3;
          fact6 = lnhx*fact4;
          fact7 = lnhy*fact5;
          fact8 = lnhz*(a1z*fact2 + a2z*fact3);
          fact9 = fact6 + fact7 + fact8;

          *hplus += fact1*(fact4*lnhx - fact5*lnhy + fact9*(e1xe1x - e1ye1y
          + e2xe2x - e2ye2y));
          *hcross += fact1*(fact4*lnhy - fact5*lnhx + fact9*(e1xe1y + e2xe2y));
        case 0:
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    case 1: // harmonic number
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
          fact1 = 1. - 2.*eta;
          fact2 = 8. + fact1*(30. + 9.*e1ze1z + 19.*e2ze2z);
          fact3 = 72. + fact1*(6. + e2ze2z - 9.*e1ze1z);
          fact4 = 40. + fact1*(18. + 15.*e2ze2z + 5.*e1ze1z);
          fact5 = 8. + fact1*(30. + 9.*e2ze2z + 19.*e1ze1z);
          fact6 = 72. + fact1*(6. + e1ze1z - 9.*e2ze2z);
          fact7 = 40. + fact1*(18. + 15.*e1ze1z + 5.*e2ze2z);
          fact8 = v5*dm/384.;

          *hplus += fact8*(((e1xe1x - e1ye1y)*e2z*fact2 - (e2xe2x
          - e2ye2y)*e2z*fact3 + 2.*e1z*(e1ye2y - e1xe2x)*fact4) + I*((-(e2xe2x
          - e2ye2y)*fact5 + (e1xe1x - e1ye1y)*fact6)*e1z - 2.*e2z*(e1ye2y
          - e1xe2x)*fact7));
          *hcross += (2.*fact8)*((-e2xe2y*e2z*fact3 + e1xe2z*e1y*fact2
          - e1z*(e1xe2y + e1ye2x)*fact4) + I*((e1xe2y + e1ye2x)*e2z*fact7
          + (e1xe1y*fact6 - e2xe2y*fact5)*e1z));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 2:
          fact1 = v4*0.25;

          *hplus += fact1*(((a2y - a1y)*e1x - (a1x - a2x)*e1y) + I*((a2y
          - a1y)*e2x - (a1x - a2x)*e2y));
          *hcross += fact1*(((a1x - a2x)*e1x - (a1y - a2y)*e1y) + I*((a1x
          - a2x)*e2x - (a1y - a2y)*e2y));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 1:
          fact1 = e1xe2x - e1ye2y;
          fact2 = e1ye1y - e1xe1x;
          fact3 = e2xe2x - e2ye2y;
          fact4 = e1xe2y + e1ye2x;
          fact5 = e1xe1y;
          fact6 = e2xe2y;
          fact7 = dm*v3*0.0625;

          *hplus += fact7*((6.*e1z*fact1 + e2z*(5.*fact2 + fact3))
          + I*(e1z*(fact2 + 5.*fact3) - 6.*e2z*fact1));
          *hcross += (2.*fact7)*((3.*e1z*fact4 + e2z*(-5.*fact5 + fact6))
          + I*(e1z*(5.*fact6 - fact5) - 3.*e2z*fact4));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 0:
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    case 2: // harmonic number
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
          logfac = log(v/v0);
          fact1 = e1xe2x - e1ye2y;
          fact2 = -e1xe1x + e1ye1y + e2xe2x - e2ye2y;
          fact3 = e1ye2x + e1xe2y;
          fact4 = -e1xe1y + e2xe2y;

          *hplus += v5*((12.*fact1*logfac + fact2*LAL_PI) + I*(6.*fact2*logfac
          - 2.*fact1*LAL_PI));
          *hcross += v5*((2.*(6.*fact3*logfac + fact4*LAL_PI))
          + I*(2.*(6.*fact4*logfac - fact3*LAL_PI)));

          fact1 = a1x*(7. + dm) + a2x*(7. - dm);
          fact2 = a1y*(7. + dm) + a2y*(7. - dm);
          fact3 = a1z*(11. - 3.*dm) + a2z*(11. + 3.*dm);
          fact4 = a1x*(41. - dm) + a2x*(41. + dm);
          fact5 = a1y*(41. - dm) + a2y*(41. + dm);
          fact6 = a1z*(41. - dm) + a2z*(41. + dm);
          fact7 = lnhx*fact4 + lnhy*fact5 + lnhz*fact6;
          fact8 = e1xe1x - e1ye1y - (e2xe2x - e2ye2y);
          fact9 = v5/48.;

          *hplus += fact9*((3.*(e1ye2z + e1ze2y)*fact1 + 3.*(e1xe2z
          + e1ze2x)*fact2 - 6.*(e1ye2x + e1xe2y)*fact3 + fact8*fact7)
          + I*(-3.*(e1ye1z - e2ye2z)*fact1 - 3.*(e1xe1z - e2xe2z)*fact2
          + 6.*(e1xe1y - e2xe2y)*fact3 + 2.*(e1xe2x - e1ye2y)*fact7));
          *hcross += fact9*((-3.*(e1ze2x + e1xe2z)*fact1 + 3.*(e1ze2y
          + e1ye2z)*fact2 + 6.*(e1xe2x - e1ye2y)*fact3 + 2.*(e1xe1y
          - e2xe2y)*fact7) + I*(3.*(e1xe1z - e2xe2z)*fact1 - 3.*(e1ye1z
          - e2ye2z)*fact2 - 3.*fact8*fact3 + 2.*(e1ye2x + e1xe2y)*fact7));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif

        case 2:
          fact5 = -1. + 3.*eta;
          fact1 = -13. + eta + (6.*e2ze2z + 2.*e1ze1z)*fact5;
          fact2 = -13. + eta + (6.*e1ze1z + 2.*e2ze2z)*fact5;
          fact3 = e1ze2z*fact5;
          fact4 = -13. + eta + 4.*(e1ze1z + e2ze2z)*fact5;
          fact6 = v4/6.;

          *hplus += fact6*((((e1ye1y - e1xe1x)*fact1 + (e2xe2x
          - e2ye2y)*fact2)*0.5) + I*(2.*(e1xe1x - e1ye1y + e2xe2x
          - e2ye2y)*fact3 + (e1ye2y - e1xe2x)*fact4));
          *hcross += fact6*((-e1xe1y*fact1 + e2xe2y*fact2) + I*(4.*(e1xe1y
          + e2xe2y)*fact3 - (e1ye2x + e1xe2y)*fact4));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
          *hplus += v2*(0.5*(e1ye1y - e2ye2y + e2xe2x - e1xe1x) + I*(e1ye2y
          - e1xe2x));
          *hcross += v2*((e2xe2y - e1xe1y) - I*(e1ye2x + e1xe2y));
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    case 3: // harmonic number
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
          fact1 = v5*dm*9./256.;
          fact2 = 1. - 2.*eta;
          fact3 = 48. + fact2*(4. + 33.*e1ze1z + 9.*e2ze2z);
          fact4 = 48. + fact2*(4. + 15.*e1ze1z + 15.*e2ze2z);
          fact5 = 48. + fact2*(4. - 3.*e1ze1z + 21.*e2ze2z);
          fact6 = 48. + fact2*(4. + 33.*e2ze2z + 9.*e1ze1z);
          fact7 = 48. + fact2*(4. - 3.*e2ze2z + 21.*e1ze1z);

          *hplus += fact1*(((e2xe2x - e2ye2y)*e2z*fact3 + 2.*e1z*(e1ye2y
          - e1xe2x)*fact4 - (e1xe1x - e1ye1y)*e2z*fact5) + I*(2.*(e1ye2y
          - e1xe2x)*e2z*fact4 + (e1xe1x - e1ye1y)*e1z*fact6 - e1z*(e2xe2x
          - e2ye2y)*fact7));
          *hcross += fact1*((2.*(e2xe2y*e2z*fact3 - (e1xe2y + e1ye2x)*e1z*fact4
          - e1xe1y*e2z*fact5)) + I*(2.*(-e1z*e2xe2y*fact7 + e1xe1y*e1z*fact6
          - (e1xe2y + e1ye2x)*e2z*fact4)));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
          fact1 = v3*dm*9./16.;
          fact2 = 2.*(e1xe2x - e1ye2y);
          fact3 = e1xe1x - e1ye1y - (e2xe2x - e2ye2y);
          fact4 = 2.*(e1xe2y + e1ye2x);
          fact5 = 2.*(e1xe1y - e2xe2y);

          *hplus += fact1*((e1z*fact2 + e2z*fact3) - I*(e1z*fact3 - e2z*fact2));
          *hcross += fact1*((e1z*fact4 + e2z*fact5) + I*(-e1z*fact5
          + e2z*fact4));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 0:
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    case 4: // harmonic number
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
        case 2:
          fact1 = v4*4.*(1. - 3.*eta)/3.;
          fact2 = e1xe2x - e1ye2y;
          fact3 = e1xe1x - e1ye1y - (e2xe2x - e2ye2y);
          fact4 = e1ze1z - e2ze2z;
          fact5 = e1xe1y - e2xe2y;
          fact6 = e1ye2x + e1xe2y;

          *hplus = fact1*((0.5*fact4*fact3 - 2.*e1ze2z*fact2) + I*(fact4*fact2
          + e1ze2z*fact3));
          *hcross = fact1*((fact4*fact5 - 2.*e1ze2z*fact6) + I*(fact4*fact6
          + 2.*e1ze2z*fact5));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    case 5: // harmonic number
      switch( ampO )
      {
        case -1: /* Use highest known PN order - move if new orders added */
        case 3:
          fact1 = -v5*dm*(1. - 2.*eta)*625./384.;
          fact2 = e1xe2x - e1ye2y;
          fact3 = e1xe1x - e1ye1y - (e2xe2x - e2ye2y);
          fact4 = e1z*(e1ze1z - 3.*e2ze2z);
          fact5 = e2z*(e2ze2z - 3.*e1ze1z);
          fact6 = e1ye2x + e1xe2y;
          fact7 = e1xe1y - e2xe2y;

          *hplus += (fact1)*((fact4*fact2 - 0.5*fact5*fact3) - I*(fact5*fact2
          + 0.5*fact4*fact3));
          *hcross += (fact1)*((fact4*fact6 - fact5*fact7) - I*(fact4*fact7
          + fact5*fact6));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
          __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d, highest is %d\n", __func__, ampO, 3 );
          break;
      }
      break;

    default: // harmonic number. zero at this order
      break;
  }
  return XLAL_SUCCESS;
}

/** @} */

/**
 * @defgroup lalsimulation_inference LALSimulation-LALInference parameter transformations
 * @author Riccardo Sturani
 *
 * @brief Functions to transform waveform parameters between LALSimulation and LALInference coordinate conventions
 */

/**
 * @ingroup lalsimulation_inference
 * @brief Transform Precessing Parameters
 *
 * @details Routine for transforming LALInference geometric variables to ChooseWaveform input
 *
 * Function to specify the desired orientation of a precessing binary in terms
 * of several angles and then compute the vector components with respect to
 * orbital angular momentum as needed to specify binary configuration for
 * ChooseTDWaveform.
 *
 * @image html lalsiminspiral_orbitelementsJ.svg Representation of input variables.
 *
 * ### Input:
 * * thetaJN is the inclination between total angular momentum (J) and the
 *   direction of propagation (N=(0,sin(thetaJN),cos(thetaJN)))
 *   @note This choice has been made to made so that thetaJN -> inclination
 *   for \f$ S_{1}+S_{2} \to 0\f$.
 * * theta1 and theta2 are the inclinations of \f$S_{1,2}\f$
 *   measured from the Newtonian orbital angular momentum (\f$L_{N}\f$).
 * * phi12 is the difference in azimuthal angles of \f$S_{1,2}\f$.
 * * chi1, chi2 are the dimensionless spin magnitudes ( chi1, chi2 \f$\le 1\f$)
 * * phiJL is the azimuthal angle of \f$L_{N}\f$ on its cone about J.
 * * m1, m2, f_ref, phiref are the component masses and reference GW frequency
     and orbital phase, they are needed to compute the magnitude of
     \f$L_{N}\f$, and thus J.
 *
 * ### Output:
 * incl - inclination angle of N relative to L_N (N=(0,sin(incl),cos(incl)))
 * in the p-q-Z frame.
 * x, y, z components \f$S_{1,2}\f$ (unit spin vectors times their
 * dimensionless spin magnitudes - i.e. they have unit magnitude for
 * extremal BHs and smaller magnitude for slower spins).
 * where x-y are rotated by phiRef with respect to p-q,
 * i.e. is \f$S_{1}\f$ wrt to x-y is (a,b,0), wrt to p-q will be
 * (a cos(phiRef) + b sin (phiRef), )
 * @note Here the \"total\" angular momentum is computed as
 * J = \f$L_{N}(1+l_{1PN}) + S_{1} + S_{2}\f$
 * where \f$L_N\f$ is the Newtonian orbital angular momentum and \f$l_{1PN}\f$
 * its relative 1PN corrections. In fact, there are
 * PN corrections to L which contribute to J that are NOT ACCOUNTED FOR
 * in this function. This is done so to avoid complications with spin-orbit
 * contributions to L, which would require the full knowledge fo the orbital
 * motion, not just the evolution of L (see e.g. eq.2.9c of
 * arXiv:gr-qc/9506022). Also, it is believed that the difference in Jhat
 * with or without these PN corrections to L is quite small.
 *
 * @attention fRef = 0 is not a valid choice. If you will pass fRef=0 into
 * ChooseWaveform, then here pass in f_min, the starting GW frequency
 *
 * UNREVIEWED
 */

int XLALSimInspiralTransformPrecessingNewInitialConditions(
        REAL8 *incl,	/**< Inclination angle of L_N (returned) */
        REAL8 *S1x,	/**< S1 x component (returned) */
        REAL8 *S1y,	/**< S1 y component (returned) */
        REAL8 *S1z,	/**< S1 z component (returned) */
        REAL8 *S2x,	/**< S2 x component (returned) */
        REAL8 *S2y,	/**< S2 y component (returned) */
        REAL8 *S2z,	/**< S2 z component (returned) */
        const REAL8 thetaJN, 	/**< zenith angle between J and N (rad) */
        const REAL8 phiJL,  	/**< azimuthal angle of L_N on its cone about J (rad) */
        const REAL8 theta1,  	/**< zenith angle between S1 and LNhat (rad) */
        const REAL8 theta2,  	/**< zenith angle between S2 and LNhat (rad) */
        const REAL8 phi12,  	/**< difference in azimuthal angle btwn S1, S2 (rad) */
        const REAL8 chi1,	/**< dimensionless spin of body 1 */
        const REAL8 chi2,	/**< dimensionless spin of body 2 */
        const REAL8 m1_SI,	/**< mass of body 1 (kg) */
        const REAL8 m2_SI,	/**< mass of body 2 (kg) */
        const REAL8 fRef,	/**< reference GW frequency (Hz) */
        const REAL8 phiRef	/**< reference orbital phase */
        )
{
    /* Check that fRef is sane */
    if( fRef == 0. )
    {
        XLALPrintError("XLAL Error - %s: fRef=0 is invalid. Please pass in the starting GW frequency instead.\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( (chi1<0.) || (chi1>1.) || (chi2<0.) || (chi2>1.) )
    {
      XLALPrintError("XLAL Error - %s: chi1,2=0  must be between 0 and 1, values %8.4f -- %8.4f passed.\n", __func__,chi1,chi2);
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 m1, m2, eta, v0, theta0, phi0, Jnorm, tmp1, tmp2;
    REAL8 Jhatx, Jhaty, Jhatz, LNhx, LNhy, LNhz, Jx, Jy, Jz, Lmag;
    REAL8 s1hatx,s1haty,s1hatz,s2hatx,s2haty,s2hatz;
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z;

    /* Starting frame: LNhat is along the z-axis and the unit
     * spin vectors are defined from the angles relative to LNhat.
     * Note that we put s1hat in the x-z plane, and phi12
     * sets the azimuthal angle of s2hat measured from the x-axis.
     */
    LNhx = 0.;
    LNhy = 0.;
    LNhz = 1.;
    /* Spins are given wrt to L,
         * but still we cannot fill the spin as we do not know
     * what will be the relative orientation of L and N.
         * Note that these spin components are NOT wrt to binary
     * separation vector, but wrt to binary separation vector
     * at phiref=0.
     */
    s1hatx = sin(theta1)*cos(phiRef);
    s1haty = sin(theta1)*sin(phiRef);
    s1hatz = cos(theta1);
    s2hatx = sin(theta2) * cos(phi12+phiRef);
    s2haty = sin(theta2) * sin(phi12+phiRef);
    s2hatz = cos(theta2);

    /* Define several internal variables needed for magnitudes */
    m1 = m1_SI/LAL_MSUN_SI;
    m2 = m2_SI/LAL_MSUN_SI;
    eta=m1*m2/(m1+m2)/(m1+m2);
    // v parameter at reference point
    v0 = cbrt( (m1+m2) * LAL_MTSUN_SI *LAL_PI * fRef );

    /* Define S1, S2, J with proper magnitudes */
    Lmag = XLALSimInspiralLN(m1+m2,eta,v0)*(1.+v0*v0*XLALSimInspiralL_2PN(eta));
    s1x = m1 * m1 * chi1 * s1hatx;
    s1y = m1 * m1 * chi1 * s1haty;
    s1z = m1 * m1 * chi1 * s1hatz;
    s2x = m2 * m2 * chi2 * s2hatx;
    s2y = m2 * m2 * chi2 * s2haty;
    s2z = m2 * m2 * chi2 * s2hatz;
    Jx = s1x + s2x;
    Jy = s1y + s2y;
    Jz = Lmag + s1z + s2z;

    /* Normalize J to Jhat, find its angles in starting frame */
    Jnorm = sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
    Jhatx = Jx / Jnorm;
    Jhaty = Jy / Jnorm;
    Jhatz = Jz / Jnorm;
    theta0 = acos(Jhatz);
    phi0 = atan2(Jhaty, Jhatx);

    /* Rotation 1: Rotate about z-axis by -phi0 to put Jhat in x-z plane */
    ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
    ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
    //do not need to perform explicitly the rotation on L and J
    //ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);
    //ROTATEZ(-phi0, LNhx, LNhy, LNhz);

    /* Rotation 2: Rotate about new y-axis by -theta0
     * to put Jhat along z-axis
     */
    ROTATEY(-theta0, LNhx, LNhy, LNhz);
    ROTATEY(-theta0, s1hatx, s1haty, s1hatz);
    ROTATEY(-theta0, s2hatx, s2haty, s2hatz);
    //do not need to perform explicitly the rotation on J
    //ROTATEY(-theta0, Jhatx, Jhaty, Jhatz);

    /* Rotation 3: Rotate about new z-axis by phiJL to put L at desired
     * azimuth about J. Note that is currently in x-z plane towards -x
     * (i.e. azimuth=pi). Hence we rotate about z by phiJL - LAL_PI
     */
    ROTATEZ(phiJL - LAL_PI, LNhx, LNhy, LNhz);
    ROTATEZ(phiJL - LAL_PI, s1hatx, s1haty, s1hatz);
    ROTATEZ(phiJL - LAL_PI, s2hatx, s2haty, s2hatz);
    //do not need to perform explicitly the rotation on J
    //ROTATEZ(phiJL - LAL_PI, Jhatx, Jhaty, Jhatz);

    /* The cosinus of the angle between L and N is the scalar
         * product of the two vectors.
         * We do not need to perform additional rotation to compute it.
         */
    REAL8 Nx=0.;
    REAL8 Ny=sin(thetaJN);
    REAL8 Nz=cos(thetaJN);
    *incl=acos(Nx*LNhx+Ny*LNhy+Nz*LNhz); //output

    /* Rotation 4-5: Now J is along z and N in y-z plane, inclined from J
         * by thetaJN and with >ve component along y.
         * Now we bring L into the z axis to get spin components.
     */
    REAL8 thetaLJ = acos(LNhz);
    REAL8 phiL    = atan2(LNhy, LNhx);

    ROTATEZ(-phiL, s1hatx, s1haty, s1hatz);
    ROTATEZ(-phiL, s2hatx, s2haty, s2hatz);
    ROTATEZ(-phiL, Nx, Ny, Nz);
    // do not need to perform explicitly the rotations on L and J
    //ROTATEZ(-phiL, LNhx, LNhy, LNhz);
    //ROTATEZ(-phiL, Jhatx, Jhaty, Jhatz);

    ROTATEY(-thetaLJ, s1hatx, s1haty, s1hatz);
    ROTATEY(-thetaLJ, s2hatx, s2haty, s2hatz);
    ROTATEY(-thetaLJ, Nx, Ny, Nz);
    // do not need to perform explicitly the rotations on L and J
    //ROTATEY(-thetaLJ, LNhx, LNhy, LNhz);
    //ROTATEY(-thetaLJ, Jhatx, Jhaty, Jhatz);

    /* Rotation 6: Now L is along z and we have to bring N
     * in the y-z plane with >ve y components.
     */
    REAL8 phiN = atan2(Ny, Nx);
    //Note the extra -phiRef here:
        // output spins must be given wrt to two body separations
        // which are rigidly rotated with spins
    ROTATEZ(LAL_PI/2.-phiN-phiRef, s1hatx, s1haty, s1hatz);
    ROTATEZ(LAL_PI/2.-phiN-phiRef, s2hatx, s2haty, s2hatz);
    // do not need to perform explicitly the rotations on L, J and N
    //ROTATEZ(LAL_PI/2.-phiN, LNhx, LNhy, LNhz);
    //ROTATEZ(LAL_PI/2.-phiN, Nx, Ny, LNz);
    //ROTATEZ(LAL_PI/2.-phiN, Jhatx, Jhaty, Jhatz);

    /* Set pointers to rotated spin vectors */
    *S1x = s1hatx*chi1;
    *S1y = s1haty*chi1;
    *S1z = s1hatz*chi1;
    *S2x = s2hatx*chi2;
    *S2y = s2haty*chi2;
    *S2z = s2hatz*chi2;

    //Uncomment the following lines for a check of the rotation
    /*printf("*****************************************\n");
    printf("** Check of TransformPrec...Conditions **\n");
    printf("*****************************************\n");
    //Rot 1:
    ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);
    //Rot 2:
    ROTATEY(-theta0, Jhatx, Jhaty, Jhatz);
    //Rot 3:
    ROTATEZ(phiJL - LAL_PI, Jhatx, Jhaty, Jhatz);
    //Rot 4:
    ROTATEZ(-phiL, Nx, Ny, Nz);
    ROTATEZ(-phiL, LNhx, LNhy, LNhz);
    ROTATEZ(-phiL, Jhatx, Jhaty, Jhatz);
    //Rot 5:
    ROTATEY(-thetaLJ, Nx, Ny, Nz);
    ROTATEY(-thetaLJ, LNhx, LNhy, LNhz);
    ROTATEY(-thetaLJ, Jhatx, Jhaty, Jhatz);
    //Rot 6:
    ROTATEZ(LAL_PI/2.-phiN, Nx, Ny, Nz);
    ROTATEZ(LAL_PI/2.-phiN, LNhx, LNhy, LNhz);
    ROTATEZ(LAL_PI/2.-phiN, Jhatx, Jhaty, Jhatz);
    printf("LNhat: %12.4e  %12.4e  %12.4e\n",LNhx,LNhy,LNhz);
    printf("       %12.4e  %12.4e  %12.4e\n",0.,0.,1.);
    printf("N:     %12.4e  %12.4e  %12.4e\n",Nx,Ny,Nz);
    printf("       %12.4e  %12.4e  %12.4e\n",0.,sin(*incl),cos(*incl));
    printf("J.Lhat  i: %12.4e f: %12.4e\n",Jz,Jnorm*Jhatz);
    printf("S1.L    i: %12.4e f: %12.4e\n",chi1*cos(theta1),*S1z);
    printf("S2.L    i: %12.4e f: %12.4e\n",chi2*cos(theta2),*S2z);
    printf("S1.S2   i: %12.4e f: %12.4e\n",chi1*chi2*(sin(theta1)*sin(theta2)*cos(phi12)+cos(theta1)*cos(theta2)),(*S1x)*(*S2x)+(*S1y)*(*S2y)+(*S1z)*(*S2z));
    printf("S1.J i: %12.4e f: %12.4e\n",chi1*(sin(theta1)*Jx+cos(theta1)*Jz),Jnorm*((*S1x)*Jhatx+(*S1y)*Jhaty+(*S1z)*Jhatz));
    printf("S2.J i: %12.4e f: %12.4e\n",chi2*(sin(theta2)*(cos(phi12)*Jx+sin(phi12)*Jy)+cos(theta2)*Jz),Jnorm*((*S2x)*Jhatx+(*S2y)*Jhaty+(*S2z)*Jhatz));
    printf("Jhat.Nhat i: %12.4e f: %12.4e\n",cos(thetaJN),Jhatx*Nx+Jhaty*Ny+Jhatz*Nz);
    printf("Norm Jhat: %12.4e  norm N: %12.4e\n",Jhatx*Jhatx+Jhaty*Jhaty+Jhatz*Jhatz,Nx*Nx+Ny*Ny+Nz*Nz);
    printf("*****************************************\n");*/
    return XLAL_SUCCESS;
}

/**
 * @ingroup lalsimulation_inference
 * @brief inverse to XLALSimInspiralTransformPrecessingNewInitialConditions()
 *
 * @details This function performs inverse transformation to
 * XLALSimInspiralTransformPrecessingNewInitialConditions()
 * it takes as input waveform parameters, assume to be defined in the
 * L=z, n=x (L-robital momentum at fRef, n is orbital separation at fRef.
 * Direction of propagation (direction to the observer) N is defined by spherical angles
 * (pi/2-phiRef, inclination).
 * The return parameters are what is used in PE for sampling (see description in ....)
 * Note that the masses are in *solar mass* and |L| is computed to the same order as in the
 * direct function above. Spins are dimensionless.
 */

/** @{ */

int XLALSimInspiralTransformPrecessingWvf2PE(
    REAL8 *thetaJN, 	/**< zenith angle between J and N (rad) [return]*/
        REAL8 *phiJL,  	/**< azimuthal angle of L_N on its cone about J (rad) [return] */
        REAL8 *theta1,  	/**< zenith angle between S1 and LNhat (rad) [return] */
        REAL8 *theta2,  	/**< zenith angle between S2 and LNhat (rad) [return] */
        REAL8 *phi12,  	/**< difference in azimuthal angle btwn S1, S2 (rad) [return] */
        REAL8 *chi1,	/**< dimensionless spin of body 1 */
        REAL8 *chi2,	/**< dimensionless spin of body 2 */
        const REAL8 incl,	/**< Inclination angle of L_N (returned) */
        const REAL8 S1x,	/**< S1 x component (input) */
        const REAL8 S1y,	/**< S1 y component (input) */
        const REAL8 S1z,	/**< S1 z component (input) */
        const REAL8 S2x,	/**< S2 x component (input) */
        const REAL8 S2y,	/**< S2 y component (input) */
        const REAL8 S2z,	/**< S2 z component (input) */
        const REAL8 m1,  	/**< mass of body 1 (solar mass) */
        const REAL8 m2,	        /**< mass of body 2 (solar mass) */
        const REAL8 fRef,	/**< reference GW frequency (Hz) */
        const REAL8 phiRef	/**< reference orbital phase */
        )
{
   /* Check that fRef is sane */
    if( fRef == 0. )
    {
        XLALPrintError("XLAL Error - %s: fRef=0 is invalid. Please pass in the starting GW frequency instead.\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 eta, v0, Jnorm, tmp1, tmp2; // theta0, phi0, Jnorm, tmp1, tmp2;
    REAL8 Jhatx, Jhaty, Jhatz, LNhx, LNhy, LNhz, Jx, Jy, Jz, Lmag;
    REAL8 s1hatx,s1haty,s1hatz,s2hatx,s2haty,s2hatz;
    REAL8 phi1, phi2, thetaJL, phiJ;
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 Nx, Ny, Nz, phiO, phiN;

    /* Starting frame: LNhat is along the z-axis and the unit
     * spin vectors are defined from the angles relative to LNhat.
     */

    LNhx = 0.;
    LNhy = 0.;
    LNhz = 1.;
    *chi1 = sqrt(S1x*S1x + S1y*S1y + S1z*S1z);
    *chi2 = sqrt(S2x*S2x + S2y*S2y + S2z*S2z);
    if ((*chi1) > 0.0){
        s1hatx = S1x/(*chi1);
        s1haty = S1y/(*chi1);
        s1hatz = S1z/(*chi1);
    }else{
        s1hatx = 0.0;
        s1haty = 0.0;
        s1hatz = 0.0;
    }

    if ((*chi2) > 0.0){
        s2hatx = S2x/(*chi2);
        s2haty = S2y/(*chi2);
        s2hatz = S2z/(*chi2);
    }else{
        s2hatx = 0.0;
        s2haty = 0.0;
        s2hatz = 0.0;
    }

    phi1 = atan2(s1haty, s1hatx);
    phi2 = atan2(s2haty, s2hatx);

    *phi12 = phi2 - phi1;
    if (*phi12 < 0.0){
        *phi12 += 2.0*LAL_PI;
    }
    *theta1 = acos(s1hatz);
    *theta2 = acos(s2hatz);

    eta=m1*m2/(m1+m2)/(m1+m2);
    // v parameter at reference point
    v0 = cbrt( (m1+m2) * LAL_MTSUN_SI *LAL_PI * fRef );

    /* Define S1, S2, J with proper magnitudes */
    Lmag = XLALSimInspiralLN(m1+m2,eta,v0)*(1.0 + v0*v0*XLALSimInspiralL_2PN(eta));
    s1x = m1 * m1 * S1x;
    s1y = m1 * m1 * S1y;
    s1z = m1 * m1 * S1z;
    s2x = m2 * m2 * S2x;
    s2y = m2 * m2 * S2y;
    s2z = m2 * m2 * S2z;
    Jx = s1x + s2x;
    Jy = s1y + s2y;
    Jz = Lmag * LNhz + s1z + s2z;

    /* Normalize J to Jhat, find its angles in starting frame */
    Jnorm = sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
    Jhatx = Jx / Jnorm;
    Jhaty = Jy / Jnorm;
    Jhatz = Jz / Jnorm;

    thetaJL = acos(Jhatz);
    phiJ = atan2(Jhaty, Jhatx);

    phiO = 0.5*LAL_PI - phiRef;
    Nx = sin(incl)*cos(phiO);
    Ny = sin(incl)*sin(phiO);
    Nz = cos(incl);

    *thetaJN = acos(Jhatx*Nx + Jhaty*Ny + Jhatz*Nz);

    /* The easiest way to define the phiJL is to rotate to the frame
     * where J is along z and N is in the y-z plane */
    ROTATEZ(-phiJ, Nx, Ny, Nz);
    ROTATEY(-thetaJL, Nx, Ny, Nz);

    ROTATEZ(-phiJ, LNhx, LNhy, LNhz);
    ROTATEY(-thetaJL, LNhx, LNhy, LNhz);
    /* You can check the rotation by uncommenting the lines below*/
    /*ROTATEZ(-phiJ, Jhatx, Jhaty, Jhatz);
    ROTATEY(-thetaJL, Jhatx, Jhaty, Jhatz);*/

    phiN = atan2(Ny, Nx);
    /* N in J-frame should be in y-z plane
     * After rotation defined below N should be in y-z plane inclined by thetaJN to J=z*/
    /*ROTATEZ(0.5*LAL_PI - phiN, Nx, Ny, Nz);*/
    ROTATEZ(0.5*LAL_PI - phiN, LNhx, LNhy, LNhz);

    *phiJL = atan2(LNhy, LNhx);

    if (*phiJL < 0.0){
        *phiJL += 2.0*LAL_PI;
    }

    /* That's all folks */
    return XLAL_SUCCESS;

}

/** @} */

/**
 * @name Routines for Handling Approximants, Order, Axis, Mode Information
 * @{
 */

/**
 * Checks whether the given approximant is implemented in lalsimulation's XLALSimInspiralChooseTDWaveform().
 *
 * returns 1 if the approximant is implemented, 0 otherwise.
 */
int XLALSimInspiralImplementedTDApproximants(
    Approximant approximant /**< post-Newtonian approximant for use in waveform production */
    )
{
    switch (approximant)
    {
        case TaylorEt:
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case TaylorT4:
	    case EccentricTD:
        case EOBNRv2:
        case HGimri:
        case IMRPhenomA:
        case EOBNRv2HM:
        case SpinTaylorT5:
        case SpinTaylorT4:
        case SpinTaylorT1:
        case IMRPhenomB:
        case PhenSpinTaylor:
        case IMRPhenomC:
        case IMRPhenomD:
        case IMRPhenomHM:
        case IMRPhenomPv2:
        case IMRPhenomPv3:
        case IMRPhenomPv3HM:
        case IMRPhenomPv2_NRTidal:
        case IMRPhenomPv2_NRTidalv2:
        case IMRPhenomNSBH:
        case IMRPhenomD_NRTidalv2:
		case IMRPhenomXAS:
        case IMRPhenomXAS_NRTidalv2:
        case IMRPhenomXAS_NRTidalv3:
		case IMRPhenomXHM:
		case IMRPhenomXP:
        case IMRPhenomXP_NRTidalv2:
        case IMRPhenomXP_NRTidalv3:
		case IMRPhenomXPHM:
        case PhenSpinTaylorRD:
        case SEOBNRv1:
        case SpinDominatedWf:
        case SEOBNRv2:
        case SEOBNRv2_opt:
        case SEOBNRv3:
        case SEOBNRv3_pert:
        case SEOBNRv3_opt:
        case SEOBNRv3_opt_rk4:
        case SEOBNRv4:
        case SEOBNRv4_opt:
        case SEOBNRv4P:
        case SEOBNRv4PHM:
        case SEOBNRv2T:
        case SEOBNRv4T:
        case SEOBNRv4_ROM_NRTidalv2_NSBH:
        case SEOBNRv4_ROM_NRTidalv2:
        case SEOBNRv5_ROM_NRTidalv3:
        case NR_hdf5:
        case NRSur7dq2:
        case NRSur7dq4:
        case TEOBResum_ROM:
        case TEOBResumS:
        case SEOBNRv4HM:
        case NRHybSur3dq8:
        case IMRPhenomT:
        case IMRPhenomTHM:
        case IMRPhenomTP:
        case IMRPhenomTPHM:
        case IMRPhenomXO4a:
        case ExternalPython:
        case SEOBNRv4HM_PA:
        case pSEOBNRv4HM_PA:
            return 1;

        default:
            return 0;
    }
}

/**
 * Checks whether the given approximant is implemented in lalsimulation's XLALSimInspiralChooseFDWaveform().
 *
 *
 * returns 1 if the approximant is implemented, 0 otherwise.
 */
int XLALSimInspiralImplementedFDApproximants(
    Approximant approximant /**< post-Newtonian approximant for use in waveform production */
    )
{
    switch (approximant)
    {
        case IMRPhenomA:
        case IMRPhenomB:
        case IMRPhenomC:
        case IMRPhenomD:
        case IMRPhenomD_NRTidal:
        case IMRPhenomD_NRTidalv2:
        case IMRPhenomNSBH:
        case IMRPhenomHM:
        case IMRPhenomP:
        case IMRPhenomPv2:
        case IMRPhenomPv2_NRTidal:
		case IMRPhenomPv2_NRTidalv2:
		case IMRPhenomXAS:
        case IMRPhenomXAS_NRTidalv2:
        case IMRPhenomXAS_NRTidalv3:
	    case IMRPhenomXHM:
		case IMRPhenomXP:
        case IMRPhenomXP_NRTidalv2:
        case IMRPhenomXP_NRTidalv3:
		case IMRPhenomXPHM:
        case EOBNRv2_ROM:
        case EOBNRv2HM_ROM:
        case SEOBNRv1_ROM_EffectiveSpin:
        case SEOBNRv1_ROM_DoubleSpin:
        case SEOBNRv2_ROM_EffectiveSpin:
        case SEOBNRv2_ROM_DoubleSpin:
        case SEOBNRv2_ROM_DoubleSpin_HI:
        case Lackey_Tidal_2013_SEOBNRv2_ROM:
        case SEOBNRv4_ROM:
        case SEOBNRv4HM_ROM:
        case SEOBNRv4_ROM_NRTidal:
        case SEOBNRv4_ROM_NRTidalv2:
        case SEOBNRv4_ROM_NRTidalv2_NSBH:
        case SEOBNRv4T_surrogate:
        case SEOBNRv5_ROM:
        case SEOBNRv5HM_ROM:
        case SEOBNRv5_ROM_NRTidalv3:
        //case TaylorR2F4:
        case TaylorF2:
        case TaylorF2Ecc:
        case TaylorF2NLTides:
        case EccentricFD:
        case SpinTaylorF2:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
        case SpinTaylorT4Fourier:
        case SpinTaylorT5Fourier:
        case NRSur4d2s:
        case IMRPhenomPv3:
        case IMRPhenomPv3HM:
        case ExternalPython:
        case IMRPhenomXO4a:
            return 1;

        default:
            return 0;
    }
}

/**
 * @brief Parses a waveform string to determine approximant, PN order, and axis choice.
 * @details
 * A waveform string can contain substrings specifying the approximant,
 * the PN order, and the frame axis.  This routine decomposes the waveform
 * string to extract this information.  Here we assume that there are no
 * extraneous characters in the waveform string that do not encode this
 * information.  If extra characters are detected then this routine returns
 * a failure code.
 *
 * If one of the output parameters is set to NULL, this routine does not
 * return the value for that parameter, and does not fail if that parameter
 * cannot be determined from the waveform string; however, the full waveform
 * string must be valid.  If the axis parameter is not NULL but information
 * about the frame axis is not found in the string then the default value
 * axis is set to the default value LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW.
 * However, if the approximant or order parameters are not NULL and
 * the approximant and order cannot be determined from the waveform string,
 * then this routine produces an error.
 *
 * Parsing is not case sensitive (using the "C" locale).
 *
 * @param[out] approximant The approximate value from Approximate enum.
 * @param[out] order The PN order value from LALPNOrder enum.
 * @param[out] axis The frame axis value from LALPNOrder enum.
 * @param[in] waveform The waveform string.
 * @retval 0 Success.
 * @retval <0 Failure.
 *
 * @note
 * Users of the SWIG-Python interface probably want to use the routines
 * XLALSimInspiralGetApproximantFromString(),
 * XLALSimInspiralGetPNOrderFromString(), and
 * XLALSimInspiralGetFrameAxisFromString()
 * since there is no way to disable required matching of the PN order
 * with the SWIG-wrapped version of this routine.
 */
int XLALSimInspiralDecomposeWaveformString(int *approximant, int *order, int *axis, const char *waveform)
{
    char *string;
    int found_approximant, found_order, found_axis;
    int failed = 0;

    if (!waveform)
        XLAL_ERROR(XLAL_EFAULT);

    string = XLALStringDuplicate(waveform);

#define DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, list) delete_substring_in_list_from_string(string, list, sizeof(list)/sizeof(*list))
    found_order       = DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, lalSimulationPNOrderNames);
    found_approximant = DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, lalSimulationApproximantNames);
    found_axis        = DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, lalSimulationFrameAxisNames);
#undef DELETE_SUBSTRING_IN_LIST_FROM_STRING

    /* assign values to output parameters */
    if (approximant) {
        *approximant = found_approximant;
        /* fail if couldn't find approximant */
        if (found_approximant < 0)
            failed = 1;
    }
    if (order) {
        *order = found_order;
        /* fail if couldn't find order */
        if (found_order < 0)
            failed = 1;
    }
    if (axis) {
        *axis = found_axis;
        /* set frame axis to view if couldn't find, but don't fail */
        if (found_axis < 0)
            *axis = LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;
    }

    /* check to see if there are extra characters */
    if (strspn(string, "\b") != strlen(string))
        failed = 1;

    XLALFree(string);

    if (failed)
        XLAL_ERROR(XLAL_EINVAL, "Invalid waveform string `%s'.", waveform);
    return 0;
}

/**
 * @brief Parses a waveform string to determine approximant.
 * @details
 * This routine uses XLALSimInspiralDecomposeWaveformString() to
 * determine the approximant from the waveform string.
 * @param[in] waveform The waveform string.
 * @return The Approximant enum value, or -1 on error.
 */
int XLALSimInspiralGetApproximantFromString(const char *waveform)
{
    int approximant = -1;
    if (XLALSimInspiralDecomposeWaveformString(&approximant, NULL, NULL, waveform)!=XLAL_SUCCESS)
    {
        approximant = -1;
        XLAL_ERROR(XLAL_EFUNC);
    }
    return approximant;
}

/**
 * @deprecated
 * Like XLALSimInspiralGetApproximantFromString() but doesn't demand that the
 * remainder of the waveform string be valid.
 */
int XLALGetApproximantFromString(const char *waveform)
{
    int approximant = -1;
    int errnum = 0;
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetApproximantFromString");
    XLAL_TRY(XLALSimInspiralDecomposeWaveformString(&approximant, NULL, NULL, waveform), errnum);
    if (errnum && errnum != XLAL_EINVAL) // pass any error other than XLAL_EINVAL
        XLAL_ERROR(errnum);
    /* fail if approximant wasn't found */
    if (approximant < 0)
        XLAL_ERROR(XLAL_EINVAL, "Cannot parse approximant from string `%s'.", waveform);
    return approximant;
}

/**
 * @brief Parses a waveform string to determine PN order.
 * @details
 * This routine uses XLALSimInspiralDecomposeWaveformString() to
 * determine the PN order from the waveform string.
 * @param[in] waveform The waveform string.
 * @return The LALPNOrder enum value, or -1 on error.
 */
int XLALSimInspiralGetPNOrderFromString(const char *waveform)
{
    int order = -1;
    if (XLALSimInspiralDecomposeWaveformString(NULL, &order, NULL, waveform) < 0)
        XLAL_ERROR(XLAL_EFUNC);
    return order;
}

/**
 * @deprecated
 * Like XLALSimInspiralGetPNOrderFromString() but doesn't demand that the
 * remainder of the waveform string be valid.
 */
int XLALGetOrderFromString(const char *waveform)
{
    int order = -1;
    int errnum = 0;
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetPNOrderFromString");
    XLAL_TRY(XLALSimInspiralDecomposeWaveformString(NULL, &order, NULL, waveform), errnum);
    if (errnum && errnum != XLAL_EINVAL) // pass any error other than XLAL_EINVAL
        XLAL_ERROR(errnum);
    /* fail if order wasn't found */
    if (order < 0)
        XLAL_ERROR(XLAL_EINVAL, "Cannot parse approximant from string `%s'.", waveform);
    return order;
}

/**
 * @brief Parses a waveform string to determine frame axis.
 * @details
 * This routine uses XLALSimInspiralDecomposeWaveformString() to
 * determine the frame axis from the waveform string.  If the
 * frame axis cannot be determined, the value
 * LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW is returned.
 * @param[in] waveform The waveform string.
 * @return The LALPNOrder enum value, or -1 on error.
 */
int XLALSimInspiralGetFrameAxisFromString(const char *waveform)
{
    int axis = -1;
    if (XLALSimInspiralDecomposeWaveformString(NULL, NULL, &axis, waveform) < 0)
        XLAL_ERROR(XLAL_EFUNC);
    return axis;
}

/**
 * @deprecated
 * Like XLALSimInspiralGetFrameAxisFromString() but doesn't demand that the
 * remainder of the waveform string be valid.
 */
int XLALGetFrameAxisFromString(const char *waveform)
{
    int axis = -1;
    int errnum = 0;
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetFrameAxisFromString");
    XLAL_TRY(XLALSimInspiralDecomposeWaveformString(NULL, NULL, &axis, waveform), errnum);
    if (errnum && errnum != XLAL_EINVAL) // pass any error other than XLAL_EINVAL
        XLAL_ERROR(errnum);
    /* if axis wasn't found, use view */
    if (axis < 0)
        axis = LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;
    return axis;
}

/**
 * @brief Parses a string to determine the LALSimInspiralApplyTaper enum value.
 * @details
 * Parses a string to determine the LALSimInspiralApplyTaper enum value.
 * Parsing is not case sensitive (using the "C" locale).
 * @param[in] string The string to be parsed.
 * @return The LALSimInspiralApplyTaper enum value, or -1 on error.
 */
int XLALSimInspiralGetTaperFromString(const char *string)
{
    const char **list = lalSimulationTaperNames;
    size_t size = sizeof(lalSimulationTaperNames)/sizeof(*lalSimulationTaperNames);
    size_t i;

    if (!string)
        XLAL_ERROR(XLAL_EFAULT);

    for (i = 0; i < size; ++i)
        if (list[i])
            if (XLALStringCaseCompare(string, list[i]) == 0) // found it
                return i;

    XLAL_ERROR(XLAL_EINVAL, "Invalid injection tapering string `%s'.", string);
}

/**
 * @deprecated
 * Use XLALSimInspiralGetTaperFromString() instead.
 */
int XLALGetTaperFromString(const char *string)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetTaperFromString");
    return XLALSimInspiralGetTaperFromString(string);
}

/**
 * @brief Parses a string to determine the LALSimInspiralModesChoice enum value.
 * @details
 * Parses a string to determine the LALSimInspiralModesChoice enum value.
 * Parsing is not case sensitive (using the "C" locale).
 * @param[in] string The string to be parsed.
 * @return The LALSimInspiralModesChoice enum value, or 0 on error.
 * @note The normal error code -1 is actually a valid mode choice
 * so this routine returns 0 (which is not a valid modes choice)
 * on error rather than -1.
 */
int XLALSimInspiralGetHigherModesFromString(const char *string)
{
    const char **list = lalSimulationModesChoiceNames;
    size_t size = sizeof(lalSimulationModesChoiceNames)/sizeof(*lalSimulationModesChoiceNames);
    size_t i;

    if (!string)
        XLAL_ERROR(XLAL_EFAULT);

    /* the "ALL" case is a special case */
    if (XLALStringCaseCompare(string, "ALL") == 0)
        return LAL_SIM_INSPIRAL_MODES_CHOICE_ALL;

    for (i = 0; i < size; ++i)
        if (list[i])
            if (XLALStringCaseCompare(string, list[i]) == 0) // found it
                return i;

    XLAL_ERROR_VAL(0, XLAL_EINVAL, "Invalid injection modes choice string `%s'.", string);
}

/**
 * @deprecated
 * Use XLALSimInspiralHigherModesFromString() instead.
 */
int XLALGetHigherModesFromString(const char *string)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetHigherModesFromString");
    return XLALSimInspiralGetHigherModesFromString(string);
}

/**
 * @brief Returns a string associated with an Approximant enum value.
 * @param[in] approximant The Approximant enum value.
 * @returns A constant string or NULL if there is an error.
 */
const char * XLALSimInspiralGetStringFromApproximant(Approximant approximant)
{
    const char *s;
    if ((int)(approximant) < 0 || (int)(approximant) >= NumApproximants)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    s = lalSimulationApproximantNames[approximant];
    if (!s)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    return s;
}

/**
 * @deprecated
 * Use XLALSimInspiralHigherModesFromString() instead.
 */
const char * XLALGetStringFromApproximant(Approximant approximant)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralGetStringFromApproximant");
    return XLALSimInspiralGetStringFromApproximant(approximant);
}

/**
 * @brief Returns a string associated with a LALPNOrder enum value.
 * @param[in] order The LALPNOrder enum value.
 * @returns A constant string or NULL if there is an error.
 */
const char * XLALSimInspiralGetStringFromPNOrder(LALPNOrder order)
{
    const char *s;
    if ((int)(order) < 0 || (int)(order) >= LAL_PNORDER_NUM_ORDER)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    s = lalSimulationPNOrderNames[order];
    if (!s)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    return s;
}

/**
 * @brief Returns a string associated with a LALSimInspiralApplyTaper enum value.
 * @param[in] taper The LALSimInspiralApplyTaper enum value.
 * @returns A constant string or NULL if there is an error.
 */
const char * XLALSimInspiralGetStringFromTaper(LALSimInspiralApplyTaper taper)
{
    const char *s;
    if ((int)(taper) < 0 || (int)(taper) >= LAL_SIM_INSPIRAL_TAPER_NUM_OPTS)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    s = lalSimulationTaperNames[taper];
    if (!s)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    return s;
}

/**
 * @brief Returns a string associated with a LALSimInspiralFrameAxis enum value.
 * @param[in] axis The LALSimInspiralFrameAxis enum value.
 * @returns A constant string or NULL if there is an error.
 */
const char * XLALSimInspiralGetStringFromFrameAxis(LALSimInspiralFrameAxis axis)
{
    const char *s;
    if ((int)(axis) < 0 || (size_t)(axis) >= sizeof(lalSimulationFrameAxisNames)/sizeof(*lalSimulationFrameAxisNames))
        XLAL_ERROR_NULL(XLAL_EINVAL);
    s = lalSimulationFrameAxisNames[axis];
    if (!s)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    return s;
}

/**
 * @brief Returns a string associated with a LALSimInspiralModesChoice enum value.
 * @param[in] modes The LALSimInspiralModesChoice enum value.
 * @returns A constant string or NULL if there is an error.
 */
const char * XLALSimInspiralGetStringFromModesChoice(LALSimInspiralModesChoice modes)
{
    const char *s;
    if (modes == LAL_SIM_INSPIRAL_MODES_CHOICE_ALL) // handle this case separately
        return "ALL";
    if ((int)(modes) < 0 || (size_t)(modes) >= sizeof(lalSimulationModesChoiceNames)/sizeof(*lalSimulationModesChoiceNames))
        XLAL_ERROR_NULL(XLAL_EINVAL);
    s = lalSimulationModesChoiceNames[modes];
    if (!s)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    return s;
}

int XLALSimInspiralGetSpinSupportFromApproximant(Approximant approx){

  SpinSupport spin_support=LAL_SIM_INSPIRAL_NUMSPINSUPPORT;
  switch (approx)
  {
    case SpinTaylor:
    case SpinTaylorFrameless:
    case SpinTaylorT1:
    case SpinTaylorT4:
    case SpinTaylorT5:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case SpinTaylorT3:
    case IMRPhenomP:
    case IMRPhenomPv2:
    case IMRPhenomPv2_NRTidal:
    case IMRPhenomPv2_NRTidalv2:
    case IMRPhenomPv3:
    case IMRPhenomPv3HM:
	case IMRPhenomXP:
    case IMRPhenomXP_NRTidalv2:
    case IMRPhenomXP_NRTidalv3:
	case IMRPhenomXPHM:
    case SpinTaylorT5Fourier:
    case SpinTaylorT4Fourier:
    case SpinDominatedWf:
    case SEOBNRv3:
    case SEOBNRv3_pert:
    case SEOBNRv3_opt:
    case SEOBNRv3_opt_rk4:
    case SEOBNRv4P:
    case SEOBNRv4PHM:
    case NR_hdf5:
    case NRSur4d2s:
    case NRSur7dq2:
    case NRSur7dq4:
    case IMRPhenomTP:
    case IMRPhenomTPHM:
    case IMRPhenomXO4a:
      spin_support=LAL_SIM_INSPIRAL_PRECESSINGSPIN;
      break;
    case SpinTaylorF2:
    case FindChirpPTF:
    case HGimri:
      spin_support=LAL_SIM_INSPIRAL_SINGLESPIN;
      break;
    case TaylorF2:
    case TaylorF2Ecc:
    case TaylorF2NLTides:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomB:
    case IMRPhenomC:
    case IMRPhenomD:
    case IMRPhenomD_NRTidal:
    case IMRPhenomD_NRTidalv2:
    case IMRPhenomNSBH:
    case IMRPhenomHM:
    case IMRPhenomXAS:
    case IMRPhenomXAS_NRTidalv2:
    case IMRPhenomXAS_NRTidalv3:
    case IMRPhenomXHM:
    case SEOBNRv1:
    case SEOBNRv2:
    case SEOBNRv4:
    case SEOBNRv2_opt:
    case SEOBNRv4_opt:
    case SEOBNRv2T:
    case SEOBNRv4T:
    case SEOBNRv4HM:
    case SEOBNRv1_ROM_EffectiveSpin:
    case SEOBNRv1_ROM_DoubleSpin:
    case SEOBNRv2_ROM_EffectiveSpin:
    case SEOBNRv2_ROM_DoubleSpin:
    case SEOBNRv2_ROM_DoubleSpin_HI:
    case Lackey_Tidal_2013_SEOBNRv2_ROM:
    case SEOBNRv4_ROM:
    case SEOBNRv4HM_ROM:
    case SEOBNRv4_ROM_NRTidal:
    case SEOBNRv4_ROM_NRTidalv2:
    case SEOBNRv4_ROM_NRTidalv2_NSBH:
    case SEOBNRv4T_surrogate:
    case SEOBNRv5_ROM:
    case SEOBNRv5HM_ROM:
    case SEOBNRv5_ROM_NRTidalv3:
    case TEOBResumS:
    case TaylorR2F4:
    case IMRPhenomFB:
    case FindChirpSP:
    case NRHybSur3dq8:
    case IMRPhenomT:
    case IMRPhenomTHM:
    case SEOBNRv4HM_PA:
    case pSEOBNRv4HM_PA:
      spin_support=LAL_SIM_INSPIRAL_ALIGNEDSPIN;
      break;
    case TaylorEt:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case EccentricTD:
    case EccentricFD:
    case IMRPhenomA:
    case EOBNRv2HM:
    case EOBNRv2HM_ROM:
    case EOBNRv2:
    case EOBNRv2_ROM:
    case EOBNR:
    case EOB:
    case IMRPhenomFA:
    case GeneratePPN:
    case TEOBResum_ROM:
      spin_support=LAL_SIM_INSPIRAL_SPINLESS;
      break;
    case ExternalPython:
      spin_support=LAL_SIM_INSPIRAL_CASEBYCASE_SPINSUPPORT;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimulation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }

    return spin_support;

}

int XLALSimInspiralGetSpinFreqFromApproximant(Approximant approx){

  SpinFreq spin_freq=LAL_SIM_INSPIRAL_NUMSPINFREQ;
  switch (approx)
  {
    case SEOBNRv3:
    case SEOBNRv3_pert:
    case SEOBNRv3_opt:
    case SEOBNRv3_opt_rk4:
    case SEOBNRv4P:
    case SEOBNRv4PHM:
      spin_freq=LAL_SIM_INSPIRAL_SPINS_FLOW;
      break;
    case SpinTaylor:
    case SpinTaylorFrameless:
    case SpinTaylorT1:
    case SpinTaylorT4:
    case SpinTaylorT5:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case SpinTaylorT3:
    case IMRPhenomP:
    case IMRPhenomPv2:
    case IMRPhenomPv3:
    case IMRPhenomPv3HM:
    case IMRPhenomPv2_NRTidal:
    case IMRPhenomPv2_NRTidalv2:
	case IMRPhenomXP:
    case IMRPhenomXP_NRTidalv2:
    case IMRPhenomXP_NRTidalv3:
	case IMRPhenomXPHM:
    case SpinTaylorT5Fourier:
    case SpinTaylorT4Fourier:
    case SpinDominatedWf:
    case NRSur4d2s:
    case NRSur7dq2:
    case NRSur7dq4:
    case SpinTaylorF2:
    case IMRPhenomTP:
    case IMRPhenomTPHM:
    case IMRPhenomXO4a:
      spin_freq=LAL_SIM_INSPIRAL_SPINS_F_REF;
      break;
    case FindChirpPTF:
    case HGimri:
    case TaylorF2:
    case TaylorF2Ecc:
    case TaylorF2NLTides:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomB:
    case IMRPhenomC:
    case IMRPhenomD:
    case IMRPhenomD_NRTidal:
    case IMRPhenomD_NRTidalv2:
    case IMRPhenomNSBH:
    case IMRPhenomHM:
    case IMRPhenomXAS:
    case IMRPhenomXAS_NRTidalv2:
    case IMRPhenomXAS_NRTidalv3:
    case IMRPhenomXHM:
    case SEOBNRv1:
    case SEOBNRv2:
    case SEOBNRv4:
    case SEOBNRv2_opt:
    case SEOBNRv4_opt:
    case SEOBNRv2T:
    case SEOBNRv4T:
    case SEOBNRv4HM:
    case SEOBNRv1_ROM_EffectiveSpin:
    case SEOBNRv1_ROM_DoubleSpin:
    case SEOBNRv2_ROM_EffectiveSpin:
    case SEOBNRv2_ROM_DoubleSpin:
    case SEOBNRv2_ROM_DoubleSpin_HI:
    case Lackey_Tidal_2013_SEOBNRv2_ROM:
    case SEOBNRv4_ROM:
    case SEOBNRv4_ROM_NRTidal:
    case SEOBNRv4_ROM_NRTidalv2:
    case SEOBNRv4_ROM_NRTidalv2_NSBH:
    case SEOBNRv4T_surrogate:
    case SEOBNRv4HM_ROM:
    case SEOBNRv5_ROM:
    case SEOBNRv5HM_ROM:
    case SEOBNRv5_ROM_NRTidalv3:
    case TaylorR2F4:
    case IMRPhenomFB:
    case FindChirpSP:
    case NRHybSur3dq8:
    case TaylorEt:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case EccentricTD:
    case EccentricFD:
    case IMRPhenomA:
    case EOBNRv2HM:
    case EOBNRv2HM_ROM:
    case EOBNRv2:
    case EOBNRv2_ROM:
    case EOBNR:
    case EOB:
    case IMRPhenomFA:
    case GeneratePPN:
    case TEOBResum_ROM:
    case IMRPhenomT:
    case IMRPhenomTHM:
	case TEOBResumS:
    case SEOBNRv4HM_PA:
    case pSEOBNRv4HM_PA:
      spin_freq=LAL_SIM_INSPIRAL_SPINS_NONPRECESSING;
      break;
    case NR_hdf5:
    case ExternalPython:
      spin_freq=LAL_SIM_INSPIRAL_SPINS_CASEBYCASE;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimulation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }

    return spin_freq;

}

int XLALSimInspiralGetAllowZeroMinFreqFromApproximant(Approximant approx){

  // Models for which LAL_SIM_INSPIRAL_ALLOW_ZERO_FMIN is set allow f_min=0,
  // which means that the full length of the waveform is returned. This means
  // that in XLALSimInspiralTD, XLALSimInspiralChooseTDWaveform is called
  // instead of XLALSimInspiralTDFromTD for these models. This also means that
  // the starting frequency is not altered (as done in XLALSimInspiralTDFromTD)
  // for these models, independent of what f_min is passed.

  AllowZeroMinFreq allow_zero_fmin = LAL_SIM_INSPIRAL_NUMZEROFMIN;
  switch (approx)
  {
    case NRSur7dq2:
    case NRSur7dq4:
      allow_zero_fmin=LAL_SIM_INSPIRAL_ALLOW_ZERO_FMIN;
      break;
    default:
      allow_zero_fmin=LAL_SIM_INSPIRAL_DISALLOW_ZERO_FMIN;
    }

    return allow_zero_fmin;

}

int XLALSimInspiralApproximantAcceptTestGRParams(Approximant approx){

  TestGRaccept testGR_accept=LAL_SIM_INSPIRAL_NUM_TESTGR_ACCEPT;
  switch (approx)
    {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF1:
    case TaylorR2F4:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case PadeT1:
    case PadeF1:
    case EOB:
    case BCV:
    case BCVSpin:
    case SpinTaylorT1:
    case SpinTaylorT5:
    case SpinTaylorT3:
    case SpinTaylorT4:
    case SpinTaylorFrameless:
    case SpinTaylor:
    case SpinQuadTaylor:
    case FindChirpSP:
    case FindChirpPTF:
    case HGimri:
    case GeneratePPN:
    case BCVC:
    case FrameFile:
    case AmpCorPPN:
    case NumRel:
    case NumRelNinja2:
    case EOBNR:
    case EOBNRv2:
    case EOBNRv2_ROM:
    case EOBNRv2HM:
    case EOBNRv2HM_ROM:
    case TEOBResum_ROM:
    case SEOBNRv1:
    case SEOBNRv2:
    case SEOBNRv2_opt:
    case SEOBNRv3:
    case SEOBNRv3_pert:
    case SEOBNRv3_opt:
    case SEOBNRv3_opt_rk4:
    case SEOBNRv4:
    case SEOBNRv4_opt:
    case SEOBNRv4P:
    case SEOBNRv4PHM:
    case SEOBNRv2T:
    case SEOBNRv4T:
    case SEOBNRv4HM:
    case SEOBNRv1_ROM_EffectiveSpin:
    case SEOBNRv1_ROM_DoubleSpin:
    case SEOBNRv2_ROM_EffectiveSpin:
    case SEOBNRv2_ROM_DoubleSpin:
    case SEOBNRv2_ROM_DoubleSpin_HI:
    case Lackey_Tidal_2013_SEOBNRv2_ROM:
    case TEOBResumS:
    case IMRPhenomA:
    case IMRPhenomB:
    case IMRPhenomFA:
    case IMRPhenomFB:
    case IMRPhenomFC:
    case IMRPhenomNSBH:
    case SpinTaylorT5Fourier:
    case SpinTaylorT4Fourier:
    case TaylorEt:
    case TaylorT4:
    case TaylorN:
    case SpinDominatedWf:
    case NR_hdf5:
    case NRSur4d2s:
    case NRSur7dq2:
    case NRSur7dq4:
    case NRHybSur3dq8:
    case IMRPhenomT:
    case IMRPhenomTHM:
    case IMRPhenomTP:
    case IMRPhenomTPHM:
    case NumApproximants:
    case SEOBNRv4HM_PA:
    case IMRPhenomXO4a:
      testGR_accept=LAL_SIM_INSPIRAL_NO_TESTGR_PARAMS;
      break;
    case TaylorF2:
    case TaylorF2Ecc:
    case TaylorF2NLTides:
    case SpinTaylorF2:
    case EccentricFD:
    case Eccentricity:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case EccentricTD:
    case SEOBNRv4_ROM:
    case SEOBNRv4HM_ROM:
    case SEOBNRv4_ROM_NRTidal:
    case SEOBNRv4_ROM_NRTidalv2:
    case SEOBNRv4_ROM_NRTidalv2_NSBH:
    case SEOBNRv4T_surrogate:
    case SEOBNRv5_ROM:
    case SEOBNRv5HM_ROM:
    case SEOBNRv5_ROM_NRTidalv3:
    case IMRPhenomC:
    case IMRPhenomD:
    case IMRPhenomP:
    case IMRPhenomPv2:
    case IMRPhenomPv2_NRTidal:
    case IMRPhenomPv2_NRTidalv2:
    case IMRPhenomD_NRTidal:
    case IMRPhenomD_NRTidalv2:
    case IMRPhenomHM:
    case IMRPhenomPv3:
    case IMRPhenomPv3HM:
    case pSEOBNRv4HM_PA:
    case IMRPhenomXAS:
    case IMRPhenomXHM:
    case IMRPhenomXP:
    case IMRPhenomXPHM:
    case IMRPhenomXAS_NRTidalv2:
    case IMRPhenomXAS_NRTidalv3:
    case IMRPhenomXP_NRTidalv2:
      testGR_accept=LAL_SIM_INSPIRAL_TESTGR_PARAMS;
      break;
    case IMRPhenomXP_NRTidalv3:
      testGR_accept=LAL_SIM_INSPIRAL_TESTGR_PARAMS;
      break;
    case ExternalPython:
			testGR_accept=LAL_SIM_INSPIRAL_CASEBYCASE_TESTGR_PARAMS;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimulation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }
  return testGR_accept;
};

/* Function for introducing Lorentz violating changes in FD phase; calculates eqns. 30 & 32 of arxiv 1110.2720 for the LV phase term in FD and multiplies to h+ and hx */
int XLALSimLorentzInvarianceViolationTerm(
                                          COMPLEX16FrequencySeries **hptilde, /**< Frequency-domain waveform h+ */
                                          COMPLEX16FrequencySeries **hctilde, /**< Frequency-domain waveform hx */
                                          REAL8 m1,                           /**< Mass 1 in solar masses */
                                          REAL8 m2,                           /**< Mass 2 in solar masses */
                                          REAL8 r,                            /**< distance in metres*/
                                          LALDict *LALparams                     /**< LAL dictionary containing accessory parameters */
                                          )
{
  REAL8 f0, f, df;
  COMPLEX16 hplus, hcross, tmpExp;
  REAL8 M, eta, zeta, dPhiPref, Mc, tmpVal;
  UINT4 len, i;
  M = m1+m2;
  eta = m1*m2/(M*M);
  Mc = M*pow(eta, 0.6);
  len = (*hptilde)->data->length;

  REAL8 lambda_eff = pow(10,XLALSimInspiralWaveformParamsLookupNonGRLIVLogLambdaEff(LALparams)); /* Effective wavelength-like parameter in phase in metres */
  REAL8 nonGR_alpha = XLALSimInspiralWaveformParamsLookupNonGRLIVAlpha(LALparams);   /* Exponent defined in terms of PN order characterising LIV*/
  REAL8 LIV_A_sign = XLALSimInspiralWaveformParamsLookupNonGRLIVASign(LALparams);   /* Sign of A determining the sign of LV phase */

  if ((*hctilde)->data->length != len) {
    XLALPrintError("Lengths of plus and cross polarization series do not agree \n");
    XLAL_ERROR(XLAL_EBADLEN);
  }

  f0 = (*hptilde)->f0;
  if ((*hctilde)->f0 != f0) {
    XLALPrintError("Starting frequencies of plus and cross polarization series do not agree \n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  df = (*hptilde)->deltaF;
  if ((*hctilde)->deltaF != df) {
    XLALPrintError("Frequency steps of plus and cross polarization series do not agree \n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  UINT4 k = 0;
  if (f0 == 0.0)
      k=1;

  if (nonGR_alpha == 1) {
    zeta = LIV_A_sign*LAL_PI*r/lambda_eff; /*Eqn. (32) of arxiv:1110.2720*/
    dPhiPref = zeta*log(LAL_PI*Mc*LAL_MTSUN_SI); /*Eqn. (31) of arxiv:1110.2720;the frequency dependence is treated below*/
    for (i=k; i<len; i++) {
      f = f0 + i*df;
      tmpExp = cexp(I*(dPhiPref + zeta*log(f)));
      hplus = (*hptilde)->data->data[i] * tmpExp;
      (*hptilde)->data->data[i] = hplus;
      hcross = (*hctilde)->data->data[i] * tmpExp;
      (*hctilde)->data->data[i] = hcross;
    }
  }
  else {
    zeta = LIV_A_sign*pow(LAL_PI, (2. - nonGR_alpha))*r*pow(Mc*LAL_MRSUN_SI, (1. - nonGR_alpha))/((1. - nonGR_alpha)*pow(lambda_eff, (2. - nonGR_alpha))); /*Eqn. (30) of arxiv:1110.2720*/
    dPhiPref = zeta*pow(LAL_PI*Mc*LAL_MTSUN_SI, (nonGR_alpha - 1.)); /*Eqn. (28) of arxiv:1110.2720;the frequency dependence is treated below*/
    for (i=k; i<len; i++) {
      f = f0 + i*df;
      tmpVal = pow(f, (nonGR_alpha - 1.));
      tmpExp=cexp(-I*dPhiPref*tmpVal);
      hplus = (*hptilde)->data->data[i] * tmpExp;
      (*hptilde)->data->data[i] = hplus;
      hcross = (*hctilde)->data->data[i] * tmpExp;
      (*hctilde)->data->data[i] = hcross;
    }

  }
  return XLAL_SUCCESS;
}

/** @} */

/**
 * @name Routines Determining Waveform Durations and Frequencies
 * @{
 */

/**
 * @brief Routine to compute an overestimate of the inspiral time from a given frequency.
 * @details
 * This routine estimates the time it will take for point-particle inspiral from a
 * specified frequency to infinite frequency.  The estimate is intended to be an
 * over-estimate, so that the true inspiral time is always smaller than the time this
 * routine returns.  To obtain this estimate, the 2pN chirp time is used where all
 * negative contributions are discarded.
 * @param fstart The starting frequency in Hz.
 * @param m1 The mass of the first component in kg.
 * @param m2 The mass of the second component in kg.
 * @param s1 The dimensionless spin of the first component.
 * @param s2 The dimensionless spin of the second component.
 * @return Upper bound on chirp time of post-Newtonian inspiral in seconds.
 */
REAL8 XLALSimInspiralChirpTimeBound(REAL8 fstart, REAL8 m1, REAL8 m2, REAL8 s1, REAL8 s2)
{
    double M = m1 + m2; // total mass
    double mu = m1 * m2 / M; // reduced mass
    double eta = mu / M; // symmetric mass ratio
    /* chi = (s1*m1 + s2*m2)/M <= max(|s1|,|s2|) */
    double chi = fabs(fabs(s1) > fabs(s2) ? s1 : s2); // over-estimate of chi
    /* note: for some reason these coefficients are named wrong...
     * "2PN" should be "1PN", "4PN" should be "2PN", etc. */
    double c0 = fabs(XLALSimInspiralTaylorT2Timing_0PNCoeff(M, eta));
    double c2 = XLALSimInspiralTaylorT2Timing_2PNCoeff(eta);
    /* the 1.5pN spin term is in TaylorT2 is 8*beta/5 [Citation ??]
     * where beta = (113/12 + (25/4)(m2/m1))*(s1*m1^2/M^2) + 2 <-> 1
     * [Cutler & Flanagan, Physical Review D 49, 2658 (1994), Eq. (3.21)]
     * which can be written as (113/12)*chi - (19/6)(s1 + s2)
     * and we drop the negative contribution */
    double c3 = (226.0/15.0) * chi;
    /* there is also a 1.5PN term with eta, but it is negative so do not include it */
    double c4 = XLALSimInspiralTaylorT2Timing_4PNCoeff(eta);
    double v = cbrt(LAL_PI * LAL_G_SI * M * fstart) / LAL_C_SI;
    return c0 * pow(v, -8) * (1.0 + (c2 + (c3 + c4 * v) * v) * v * v);
}

/**
 * @brief Routine to compute an overestimate of the merger time.
 * @details
 * This routine provides an upper bound on the time it will take for compact
 * binaries to plunge and merge at the end of the quasi-stationary inspiral.
 * This is quite vague since the notion of a innermost stable circular orbit
 * is ill-defined except in a test mass limit. Nevertheless, this routine
 * assumes (i) that the innermost stable circular orbit occurs at v = c / 3,
 * or r = 9 G M / c^3 (in Boyer-Lindquist coordinates), which is roughly right
 * for an extreme Kerr black hole counter-rotating with a test particle,
 * and (ii) the plunge lasts for a shorter period than one cycle at this
 * innermost stable circular orbit.
 * @param m1 The mass of the first component in kg.
 * @param m2 The mass of the second component in kg.
 * @return Upper bound on the merger time in seconds.
 */
REAL8 XLALSimInspiralMergeTimeBound(REAL8 m1, REAL8 m2)
{
    const double norbits = 1;
    double M = m1 + m2; // total mass
    double r = 9.0 * M * LAL_MRSUN_SI / LAL_MSUN_SI;
    double v = LAL_C_SI / 3.0;
    return norbits * (2.0 * LAL_PI * r / v);
}

/**
 * @brief Routine to compute an overestimate of the ringdown time.
 * @details
 * This routine provides an upper bound on the time it will take for a
 * black hole produced by compact binary merger to ring down through
 * quasinormal mode radiation.  An approximate formula for the frequency
 * and quality of the longest-lived fundamental (n=0) dominant (l=m=2)
 * quasinormal mode * is given by Eqs. (E1) and (E2) along with Table VIII
 * of Berti, Cardoso, and Will, Physical Review D (2006) 064030.
 * Waveform generators produce 10 e-folds of ringdown radiation, so
 * this routine goes up to 11 to provide an over-estimate of the ringdown time.
 * @param M The mass of the final black hole in kg.
 * @param s The dimensionless spin of the final black hole.
 * @return Upper bound on the merger time in seconds.
 * @see Emanuele Berti, Vitor Cardoso, and Clifford M. Will, Physical Review D 73, 064030 (2006) DOI: 10.1103/PhysRevD.73.064030
 */
REAL8 XLALSimInspiralRingdownTimeBound(REAL8 M, REAL8 s)
{
    const double nefolds = 11; /* waveform generators only go up to 10 */

    /* these values come from Table VIII of Berti, Cardoso, and Will with n=0, m=2 */
    const double f1 = +1.5251;
    const double f2 = -1.1568;
    const double f3 = +0.1292;
    const double q1 = +0.7000;
    const double q2 = +1.4187;
    const double q3 = -0.4990;

    double omega = (f1 + f2 * pow(1.0 - s, f3)) / (M * LAL_MTSUN_SI / LAL_MSUN_SI);
    double Q = q1 + q2 * pow(1.0 - s, q3);
    double tau = 2.0 * Q / omega; /* see Eq. (2.1) of Berti, Cardoso, and Will */
    return nefolds * tau;
}

/**
 * @brief Routine to compute an overestimate of a final black hole dimensionless spin.
 * @details
 * This routine provides an upper bound on the dimensionless spin of a black
 * hole produced in a compact binary merger.  Uses the formula in Tichy and
 * Marronetti, Physical Review D 78 081501 (2008), Eq. (1) and Table 1, for
 * equal mass black holes, or the larger of the two spins (which covers the
 * extreme mass case).  If the result is larger than a maximum realistic
 * black hole spin, truncate at this maximum value.
 * @param S1z The z-component of the dimensionless spin of body 1.
 * @param S2z The z-component of the dimensionless spin of body 2.
 * @return Upper bound on final black hole dimensionless spin.
 * @see Tichy and Marronetti, Physical Review D 78 081501 (2008).
 * @todo It has been suggested that Barausse, Rezzolla (arXiv: 0904.2577) is
 * more accurate
 */
REAL8 XLALSimInspiralFinalBlackHoleSpinBound(REAL8 S1z, REAL8 S2z)
{
    const double maximum_black_hole_spin = 0.998;
    double s;
    /* lower bound on the final plunge, merger, and ringdown time here the
     * final black hole spin is overestimated by using the formula in Tichy and
     * Marronetti, Physical Review D 78 081501 (2008), Eq. (1) and Table 1, for
     * equal mass black holes, or the larger of the two spins (which covers the
     * extreme mass case) */
    /* TODO: it has been suggested that Barausse, Rezzolla (arXiv: 0904.2577)
     * is more accurate */
    s = 0.686 + 0.15 * (S1z + S2z);
    if (s < fabs(S1z))
        s = fabs(S1z);
    if (s < fabs(S2z))
        s = fabs(S2z);
    /* it is possible that |S1z| or |S2z| >= 1, but s must be less than 1
     * (0th law of thermodynamics) so provide a maximum value for s */
    if (s > maximum_black_hole_spin)
        s = maximum_black_hole_spin;
    return s;
}

/**
 * @brief Routine to compute an underestimate of the starting frequency for a given chirp time.
 * @details
 * This routine estimates a start frequency for a binary inspiral from which the
 * actual inspiral chirp time will be shorter than the specified chirp time.
 * To obtain this estimate, only the leading order Newtonian coefficient is used.
 * The frequency returned by this routine is guaranteed to be less than the frequency
 * passed to XLALSimInspiralChirpTimeBound() if the returned value of that routine
 * is passed to this routine as tchirp.
 * @param tchirp The chirp time of post-Newtonian inspiral s.
 * @param m1 The mass of the first component in kg.
 * @param m2 The mass of the second component in kg.
 * @return Lower bound on the starting frequency of a post-Newtonian inspiral in Hz.
 */
REAL8 XLALSimInspiralChirpStartFrequencyBound(REAL8 tchirp, REAL8 m1, REAL8 m2)
{
    double M = m1 + m2; // total mass
    double mu = m1 * m2 / M; // reduced mass
    double eta = mu / M; // symmetric mass ratio
    double c0 = XLALSimInspiralTaylorT3Frequency_0PNCoeff(M);
    return c0 * pow(5.0 * M * (LAL_MTSUN_SI / LAL_MSUN_SI) / (eta * tchirp), 3.0 / 8.0);
}

/**
 * Function that gives the value of the desired frequency given some physical parameters
 *
 */
double XLALSimInspiralGetFrequency(
    REAL8 m1,                   /**< mass of companion 1 (kg) */
    REAL8 m2,                   /**< mass of companion 2 (kg) */
    const REAL8 S1x,            /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,            /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,            /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,            /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,            /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,            /**< z-component of the dimensionless spin of object 2 */
    FrequencyFunction freqFunc  /**< name of the function to use */
    )
{

    double freq; /* The return value */

    /* Variables needed for fIMRPhenom(x) */
    double chi;
    double m1Msun = m1 / LAL_MSUN_SI;
    double m2Msun = m2 / LAL_MSUN_SI;

    /* Variables needed for f(S)EOBNRv(x) */
    UINT4 SpinAlignedEOBVersion;
    Approximant approximant;
    REAL8 spin1[3];
    REAL8 spin2[3];
    int     modeL;
    int     modeM;
    COMPLEX16Vector modefreqVec;
    COMPLEX16      modeFreq;

    switch (freqFunc)
    {
        case fSchwarzISCO:
            /* Schwarzschild ISCO */
            freq = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);
            break;
        case fIMRPhenomAFinal:
            freq = XLALSimIMRPhenomAGetFinalFreq(m1Msun, m2Msun);
            break;
        case fIMRPhenomBFinal:
            chi = XLALSimIMRPhenomBComputeChi(m1Msun, m2Msun, S1z, S2z);
            freq = XLALSimIMRPhenomBGetFinalFreq(m1Msun, m2Msun, chi);
            break;
        case fIMRPhenomCFinal:
            chi = XLALSimIMRPhenomBComputeChi(m1Msun, m2Msun, S1z, S2z);
            freq = XLALSimIMRPhenomCGetFinalFreq(m1Msun, m2Msun, chi);
            break;
        case fIMRPhenomDPeak:
            freq = XLALIMRPhenomDGetPeakFreq(m1Msun, m2Msun, S1z, S2z);
            break;

        /* EOBNR ringdown frequencies all come from the same code,
         * just with different inputs */
        case fEOBNRv2HMRD:
        case fEOBNRv2RD:
        case fSEOBNRv1RD:
        case fSEOBNRv2RD:
        case fSEOBNRv4RD:
            // FIXME: Probably shouldn't hard code the modes.
            if ( freqFunc == fEOBNRv2HMRD )
            {
                modeL = 5;
                modeM = 5;
                approximant = EOBNRv2HM;
            }
            else
            {
                modeL = 2;
                modeM = 2;
                if (freqFunc == fEOBNRv2RD) approximant = EOBNRv2;
                if (freqFunc == fSEOBNRv1RD) approximant = SEOBNRv1;
                if (freqFunc == fSEOBNRv2RD) approximant = SEOBNRv2;
                if (freqFunc == fSEOBNRv4RD) approximant = SEOBNRv4;
            }
            if ( freqFunc == fEOBNRv2RD || freqFunc == fEOBNRv2HMRD )
            {
                /* Check that spins are zero */
                if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                {
                    XLALPrintError("Non-zero spins were given, but EOBNRv2 ringdown frequencies do not depend on spin.\n");
                    XLAL_ERROR(XLAL_EINVAL);
                spin1[0] = 0.; spin1[1] = 0.; spin1[2] = 0.;
                spin2[0] = 0.; spin2[1] = 0.; spin2[2] = 0.;
                }
            }
            else
            {
                spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
                spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            }

            modefreqVec.length = 1;
            modefreqVec.data   = &modeFreq;
            if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1Msun, m2Msun, spin1, spin2, modeL, modeM, 1, approximant) != XLAL_SUCCESS )
            {
                XLAL_ERROR( XLAL_EFUNC );
            }

            freq = creal(modeFreq) / (2 * LAL_PI);
            break;

        case fSEOBNRv5RD:
            modeL = 2;
            modeM = 2;
            approximant = SEOBNRv5_ROM;
            spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
            spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            modefreqVec.length = 1;
            modefreqVec.data   = &modeFreq;
            if ( XLALSimIMREOBGenerateQNMFreqV5( &modefreqVec, m1Msun, m2Msun, spin1, spin2, modeL, modeM, 1, approximant) != XLAL_SUCCESS )
            {
                XLAL_ERROR( XLAL_EFUNC );
            }

            freq = creal(modeFreq) / (2 * LAL_PI);
            break;

        case fSEOBNRv1Peak:
        case fSEOBNRv2Peak:
        case fSEOBNRv4Peak:
        case fSEOBNRv5Peak:
            if ( freqFunc == fSEOBNRv1Peak ) SpinAlignedEOBVersion = 1;
            if ( freqFunc == fSEOBNRv2Peak ) SpinAlignedEOBVersion = 2;
            if ( freqFunc == fSEOBNRv4Peak ) SpinAlignedEOBVersion = 4;
            if ( freqFunc == fSEOBNRv5Peak ) SpinAlignedEOBVersion = 5;
            freq = XLALSimIMRSpinAlignedEOBPeakFrequency(m1, m2, S1z, S2z,
                    SpinAlignedEOBVersion);
            break;
        case fTEOBResumSFinal: // MA: Replace with TEOB-related RD frequency!
                               // CAUTION: different function for BNS/NSBH/BBH cases?
            fprintf(stdout, "Final frequency for TEOBResumS not implemented yet.\n");
            freq = XLALSimIMRSpinAlignedEOBPeakFrequency(m1, m2, S1z, S2z, 2);
            break;
        default:
            XLALPrintError("Unsupported approximant\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    return freq;
}


/**
 * Function that gives the default ending frequencies of the given approximant.
 *
 */
double XLALSimInspiralGetFinalFreq(
    REAL8 m1,                               /**< mass of companion 1 (kg) */
    REAL8 m2,                               /**< mass of companion 2 (kg) */
    const REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
    )
{
    FrequencyFunction freqFunc;

    /* input conditions */
    switch (approximant)
    {
        case EccentricTD:
        case EccentricFD:
        case EOBNRv2HM:
        case EOBNRv2:
        case IMRPhenomA:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            break;

        case TaylorF2RedSpinTidal:
        case SEOBNRv1:
        case SEOBNRv2:
        case SEOBNRv2_opt:
        case SEOBNRv4:
        case SEOBNRv4_opt:
        case IMRPhenomB:
        case IMRPhenomC:
        case TEOBResumS:
            /* Check that the transverse spins are zero */
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
                XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            break;

        default:
            break;
    }

    /* select the frequency function that is associated with each approximant */
    switch (approximant)
    {
        /* non-spinning inspiral-only models */
        // CHECKME: do they really all use Schwarzschild ISCO? */
        case TaylorEt:
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case TaylorT4:
        case EccentricTD:
        case EccentricFD:
        case TaylorF2:
        case TaylorF2Ecc:
        case TaylorF2NLTides:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
            /* Schwarzschild ISCO */
            freqFunc = fSchwarzISCO;
            break;

        /* IMR models */
        case EOBNRv2HM:
            freqFunc = fEOBNRv2HMRD;
            break;

        case EOBNRv2:
            freqFunc = fEOBNRv2RD;
            break;

        case SEOBNRv1:
            freqFunc = fSEOBNRv1RD;
            break;

        case SEOBNRv2:
        case SEOBNRv2_opt:
            freqFunc = fSEOBNRv2RD;
            break;

        case SEOBNRv4:
        case SEOBNRv4_opt:
            freqFunc = fSEOBNRv4RD;
            break;

        case SEOBNRv5_ROM:
            freqFunc = fSEOBNRv5RD;
            break;

        case IMRPhenomA:
            freqFunc = fIMRPhenomAFinal;
            break;

        case IMRPhenomB:
            freqFunc = fIMRPhenomBFinal;
            break;

        case IMRPhenomC:
            freqFunc = fIMRPhenomCFinal;
            break;

        case TEOBResumS:
           freqFunc = fTEOBResumSFinal;
       break;

        // FIXME: Following I don't know how to calculate */
        /* Spinning inspiral-only time domain */
        case SpinTaylorT5:
        case SpinTaylorT4:
        case SpinTaylorT1:
        case PhenSpinTaylor:
        /* Spinning with ringdown attachment */
        case PhenSpinTaylorRD:
        /* Spinning inspiral-only frequency domain */
        case SpinTaylorF2:
        /* NR waveforms */
        case NR_hdf5:
        case NRSur4d2s:
            XLALPrintError("I don't know how to calculate final freq. for this approximant, sorry!\n");
            XLAL_ERROR(XLAL_EINVAL);
            break;

        default:
            XLALPrintError("Unsupported approximant\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    return XLALSimInspiralGetFrequency(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, freqFunc);
}

/** @} */

/**
 * @name Waveform Conditioning Helper Routines
 * @{
 */

/**
 * @brief First stage of conditioning of time-domain waveforms.
 * @details
 * The conditioning of time-domain waveforms is done in two stages:
 *
 * 1. Take a waveform that is generated so that begins a time at least textra
 * before it reaches f_min and apply a taper over this duration textra; follow
 * this up by high-pass filtering the data at frequency f_min; finally, if the
 * original waveform was zero-padded, strip off this padding.
 *
 * 2. The high-pass filtered waveform might have transients remaining at the
 * beginning and at the end; furthermore, the waveform might end at a non-zero
 * value (especially for non-IMR waveforms).  The second stage tapers one
 * cycle at f_min from the beginning of the waveform and one cycle at f_max
 * from the end of the waveform.  There is a minimum number of samples that
 * will be used in the tapering and if the waveform is shorter than twice
 * this minimum number then no Stage 2 conditioning is done.
 *
 * This routine performs Stage 1 conditioning.  This stage is only performed
 * for waveforms originally produced in the time domain.  (Frequency domain
 * waveforms that have been transformed into the time domain have a different
 * Stage 1 conditioning.)
 *
 * @param hplus  Pointer to the plus-polarization timeseries.
 * @param hcross Pointer to the cross-polarization timeseries.
 * @param textra Extra time at beginning of waveform to taper (s).
 * @param f_min  Minimum frequency for high-pass filtering (Hz).
 * @retval  0 Success.
 * @retval <0 Failure.
 */
int XLALSimInspiralTDConditionStage1(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, REAL8 textra, REAL8 f_min)
{
    size_t nzeros;
    size_t ntaper;
    size_t j;

    /* some generators zero-pad the end of the waveform: will remove this */
    nzeros = 0;
    while (hplus->data->data[hplus->data->length - nzeros - 1] == 0.0 && hcross->data->data[hcross->data->length - nzeros - 1] == 0.0)
        ++nzeros;

    /* apply tapers over the extra duration at the beginning */
    ntaper = round(textra / hplus->deltaT);
    for (j = 0; j < ntaper; ++j) {
        double w = 0.5 - 0.5 * cos(j * LAL_PI / ntaper);
        hplus->data->data[j] *= w;
        hcross->data->data[j] *= w;
    }

    /* apply time domain filter at f_min */
    XLALHighPassREAL8TimeSeries(hplus, f_min, 0.99, 8);
    XLALHighPassREAL8TimeSeries(hcross, f_min, 0.99, 8);

    /* now take off the zero padded end */
    if (nzeros) {
        XLALShrinkREAL8TimeSeries(hplus, 0, hplus->data->length - nzeros);
        XLALShrinkREAL8TimeSeries(hcross, 0, hcross->data->length - nzeros);
    }

    return 0;
}

/**
 * @brief Second stage of conditioning of time-domain waveforms.
 * @details
 * The conditioning of time-domain waveforms is done in two stages:
 *
 * 1. Take a waveform that is generated so that begins a time at least textra
 * before it reaches f_min and apply a taper over this duration textra; follow
 * this up by high-pass filtering the data at frequency f_min; finally, if the
 * original waveform was zero-padded, strip off this padding.
 *
 * 2. The high-pass filtered waveform might have transients remaining at the
 * beginning and at the end; furthermore, the waveform might end at a non-zero
 * value (especially for non-IMR waveforms).  The second stage tapers one
 * cycle at f_min from the beginning of the waveform and one cycle at f_max
 * from the end of the waveform.  There is a minimum number of samples that
 * will be used in the tapering and if the waveform is shorter than twice
 * this minimum number then no Stage 2 conditioning is done.
 *
 * This routine performs Stage 2 conditioning.  This stage is performed both
 * for waveforms originally produced in the time domain, and for waveforms that
 * were originally produced in the frequency domain and then transformed to
 * the time domain. This stage follows some form of Stage 1 conditioning,
 * which is different depending on whether the original waveform was generated
 * in the time domain or in the frequency domain.
 *
 * @param hplus  Pointer to the plus-polarization timeseries.
 * @param hcross Pointer to the cross-polarization timeseries.
 * @param f_min  Minimum frequency for tapering at the start (Hz).
 * @param f_max  Minimum frequency for tapering at the end (Hz).
 * @retval  0 Success.
 * @retval <0 Failure.
 */
int XLALSimInspiralTDConditionStage2(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, REAL8 f_min, REAL8 f_max)
{
    const size_t min_taper_samples = 4;
    size_t ntaper;
    size_t j;

    /* final tapering at the beginning and at the end */
    /* if this waveform is shorter than 2*min_taper_samples, do nothing */
    if (hplus->data->length < 2 * min_taper_samples) {
        XLAL_PRINT_WARNING("waveform is too shorter than %zu samples: no final tapering applied", 2 * min_taper_samples);
        return 0;
    }

    /* taper end of waveform: 1 cycle at f_max; at least min_taper_samples
     * note: this tapering is done so the waveform goes to zero at the next
     * point beyond the end of the data */
    ntaper = round(1.0 / (f_max * hplus->deltaT));
    if (ntaper < min_taper_samples)
        ntaper = min_taper_samples;
    for (j = 1; j < ntaper; ++j) {
        double w = 0.5 - 0.5 * cos(j * LAL_PI / ntaper);
        hplus->data->data[hplus->data->length - j] *= w;
        hcross->data->data[hcross->data->length - j] *= w;
    }

    /* there could be a filter transient at the beginning too: we should have
     * some safety there owing to the fact that we are starting at a lower
     * frequency than is needed, so taper off one cycle at the low frequency */
    ntaper = round(1.0 / (f_min * hplus->deltaT));
    if (ntaper < min_taper_samples)
        ntaper = min_taper_samples;
    for (j = 0; j < ntaper; ++j) {
        double w = 0.5 - 0.5 * cos(j * LAL_PI / ntaper);
        hplus->data->data[j] *= w;
        hcross->data->data[j] *= w;
    }

    return 0;
}


/**
 * @brief Function for determining the starting frequency
 * of the (2,2) mode when the highest order contribution starts at fLow.
 * @details
 * Compute the minimum frequency for waveform generation
 *  using amplitude orders above Newtonian.  The waveform
 *  generator turns on all orders at the orbital
 *  associated with fMin, so information from higher
 *  orders is not included at fLow unless fMin is
 *  sufficiently low.
 *
 * @param fLow  Requested lower frequency.
 * @param ampOrder Requested amplitude order.
 * @param approximant LALApproximant
 * @retval fStart The lower frequency to use to include corrections.
 */
REAL8 XLALSimInspiralfLow2fStart(REAL8 fLow, INT4 ampOrder, INT4 approximant)
{
  if (ampOrder == -1) {
      if (approximant == SpinTaylorT5 || approximant == SpinTaylorT4)
          ampOrder = MAX_PRECESSING_AMP_PN_ORDER;
      else
          ampOrder = MAX_NONPRECESSING_AMP_PN_ORDER;
  }

    REAL8 fStart;
    fStart = fLow * 2./(ampOrder+2);
    return fStart;
}

/**
 * @deprecated Use XLALSimInspiralChooseTDWaveform() instead
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the time domain.
 *
 * The parameters passed must be in SI units.
 */
int XLALSimInspiralChooseTDWaveformOLD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    const REAL8 m1,                             /**< mass of companion 1 (kg) */
    const REAL8 m2,                             /**< mass of companion 2 (kg) */
    const REAL8 S1x,                            /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                            /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                            /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                            /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                            /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                            /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                       /**< distance of source (m) */
    const REAL8 inclination,                    /**< inclination of source (rad) */
    const REAL8 phiRef,                         /**< reference orbital phase (rad) */
    const REAL8 longAscNodes,                   /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 eccentricity,                   /**< eccentrocity at reference epoch */
    const REAL8 UNUSED meanPerAno,              /**< mean anomaly of periastron */
    // frequency sampling parameters, no default value
    const REAL8 deltaT,                         /**< sampling interval (s) */
    const REAL8 f_min,                          /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    const REAL8 lambda1,                        /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    const REAL8 lambda2,                        /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    const REAL8 dQuadParam1,                    /**< (quad-monop parameter of mass 1) / m1^5 -1 (dimensionless), give 0 for BHs */
    const REAL8 dQuadParam2,                    /**< (quad-monop parameter of mass 2) / m2^5 -1 (dimensionless), give 0 for BHs */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    const int phaseO,                           /**< twice post-Newtonian order */
    const Approximant approximant               /**< post-Newtonian approximant to use for waveform production */
    )
{
    LALDict *LALparams = NULL;
    REAL8 LNhatx, LNhaty, LNhatz, E1x, E1y, E1z;
    char *numrel_data_path;
    //REAL8 tmp1, tmp2;
    int ret;
    /* N.B. the quadrupole of a spinning compact body labeled by A is
     * Q_A = - quadparam_A chi_A^2 m_A^3 (see gr-qc/9709032)
     * where quadparam = 1 for BH ~= 4-8 for NS.
     * This affects the quadrupole-monopole interaction.
     */
    REAL8 v0 = 1., quadparam1 = 1.+dQuadParam1, quadparam2 = 1.+dQuadParam2;

    /* General sanity checks that will abort
     *
     * If non-GR approximants are added, include them in
     * XLALSimInspiralApproximantAcceptTestGRParams()
     */
    if( nonGRparams && XLALSimInspiralApproximantAcceptTestGRParams(approximant) != LAL_SIM_INSPIRAL_TESTGR_PARAMS )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Support variables for precessing wfs*/
    REAL8 incl;

    UINT4 PrecEOBversion;

    /* SEOBNR flag for model version. 1 for SEOBNRv1, 2 for SEOBNRv2 */
    UINT4 SpinAlignedEOBversion;
    REAL8 spin1x,spin1y,spin1z;
    REAL8 spin2x,spin2y,spin2z;
    REAL8 polariz=longAscNodes;
    //LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if( deltaT < 1./16385. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( S1x*S1x + S1y*S1y + S1z*S1z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if( S2x*S2x + S2y*S2y + S2z*S2z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    switch (approximant)
    {
        /* non-spinning inspiral-only models */
        case TaylorEt:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, distance, inclination, amplitudeO, phaseO);
            break;

        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2,
                    0, amplitudeO, phaseO);
            break;

        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT2PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2,
                    0, amplitudeO, phaseO);
            break;

        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2,
                    0, amplitudeO, phaseO);
            break;

        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2,
                    0, amplitudeO, phaseO);
            break;

    case EccentricTD:
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
        XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralEccentricTDPNGenerator(hplus, hcross, phiRef,
            deltaT, m1, m2, f_min, f_ref, distance, inclination, eccentricity,
            amplitudeO, phaseO);
        if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
        break;

        /* non-spinning inspiral-merger-ringdown models */
        case IMRPhenomA:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomAGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, 0., distance, inclination);
            break;

        case EOBNRv2HM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // FIXME: need to create a function to take in different modes or produce an error if all modes not given
            ret = XLALSimIMREOBNRv2AllModes(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, distance, inclination);
            break;

        case EOBNRv2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMREOBNRv2DominantMode(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, distance, inclination);
            break;

        /* spinning inspiral-only models */
        case SpinTaylorT5:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(incl);
            LNhaty = 0.;
            LNhatz = cos(incl);
            E1x = 0.;
            E1y = 1.;
            E1z = 0.;
        polariz+=LAL_PI/2.;
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT5(hplus, hcross, phiRef, deltaT,
                          m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                          LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, NULL);
            break;

        // need to make a consistent choice for SpinTaylorT4 and PSpinInspiralRD waveform inputs
        // proposal: TotalJ frame of PSpinInspiralRD
        // inclination denotes the angle between the view direction
        // and J (J is constant during the evolution, J//z, both N and initial
        // L are in the x-z plane) and the spin coordinates are given wrt
        // initial ** L **.
        case SpinTaylorT4:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
            XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(incl);
            LNhaty = 0.;
            LNhatz = cos(incl);
            E1x = 0.;
            E1y = 1.;
            E1z = 0.;
        polariz+=LAL_PI/2.;
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT4OLD(hplus, hcross, phiRef, 1.,deltaT,
                    m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
            quadparam1, quadparam2, NULL,
                    phaseO, amplitudeO);
            break;

     case SpinTaylorT1:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(incl);
            LNhaty = 0.;
            LNhatz = cos(incl);
            E1x = 0.;
            E1y = 1.;
            E1z = 0.;
        polariz+=LAL_PI/2.;
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT1OLD(hplus, hcross, phiRef, 1., deltaT,
                          m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2, NULL,
                    phaseO, amplitudeO);
            break;

        case SpinDominatedWf:
                // waveform specific sanity checks
                if (S2x != 0. || S2y != 0. || S2z != 0.){
                XLALPrintError("XLAL Error : The spindominatedwf approximant is only for 1 spin case.\n");
                XLAL_ERROR(XLAL_EDOM);
                }
                /*Maximal PN amplitude order is 1.5, maximal phase order is 2 PN*/
                if (amplitudeO > 3) {
                XLALPrintError("XLAL Error : Foe the spindominatedwf approximant maximal amplitude correction is 1.5 PN\n");
                XLAL_ERROR(XLAL_EDOM);
                }
                if (phaseO > 4){
                XLALPrintError("XLAL Error : For the spindominatedwf approximant maximal phase correction is 2 PN\n");
                XLAL_ERROR(XLAL_EDOM);
                }
                incl=inclination;
                LNhatx = 0.;
                LNhaty = 0.;
                LNhatz = 1.;
                /* Call the waveform driver routine */
                ret = XLALSimInspiralSpinDominatedWaveformInterfaceTD(hplus, hcross, deltaT, m1, m2, f_min, f_ref, distance, S1x, S1y, S1z, LNhatx, LNhaty, LNhatz, incl, phaseO, amplitudeO, phiRef);
                break;

        /* spin aligned inspiral-merger-ringdown models */
        case IMRPhenomB:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomBGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, 0., distance, inclination);
            break;

        case PhenSpinTaylor:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
        polariz+=LAL_PI/2.;
            ret = XLALSimSpinInspiralGenerator(hplus, hcross, phiRef,
                           deltaT, m1, m2, f_min, f_ref, distance, incl, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                           phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, NULL);
        break;

        case IMRPhenomC:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomCGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, 0., distance, inclination, NULL);
            break;

    case IMRPhenomD:
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        // generate TD waveforms with zero inclincation so that amplitude can be
        // calculated from hplus and hcross, apply inclination-dependent factors
        // in loop below
            /* FIXME: BUSTED -- EXTRA PARAMS NOT IMPLEMENTED */
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, 0.0, phiRef, 0.0, 0.0, 0.0, deltaT, f_min, f_ref, NULL, approximant);
        REAL8 maxamp=0;
        REAL8TimeSeries *hp = *hplus;
        REAL8TimeSeries *hc = *hcross;
        INT4 maxind=hp->data->length - 1;
        INT4 loopi;
        const REAL8 cfac=cos(inclination);
        const REAL8 pfac = 0.5 * (1. + cfac*cfac);

        for (loopi=hp->data->length - 1; loopi > -1; loopi--)
        {
            REAL8 ampsqr = (hp->data->data[loopi])*(hp->data->data[loopi]) +
               (hc->data->data[loopi])*(hc->data->data[loopi]);
            if (ampsqr > maxamp)
            {
                maxind=loopi;
                maxamp=ampsqr;
            }
            hp->data->data[loopi] *= pfac;
            hc->data->data[loopi] *= cfac;
        }
        XLALGPSSetREAL8(&(hp->epoch), (-1.) * deltaT * maxind);
        XLALGPSSetREAL8(&(hc->epoch), (-1.) * deltaT * maxind);
        break;

    case IMRPhenomPv2:
            /* FIXME: BUSTED -- EXTRA PARAMS NOT IMPLEMENTED */
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, 0.0, 0.0, 0.0, deltaT, f_min, f_ref, NULL, approximant);
        break;

        case PhenSpinTaylorRD:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at the start.\n", __func__);
            /* Call the waveform driver routine */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
        polariz+=LAL_PI/2.;
            ret = XLALSimIMRPhenSpinInspiralRDGenerator(hplus, hcross, phiRef,
                            deltaT, m1, m2, f_min, f_ref, distance, incl, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                            phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, NULL);
            break;

        case SEOBNRv1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 1;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, LALparams);
            break;

        case SEOBNRv2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 2;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, LALparams);
            break;

        case SEOBNRv4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 4;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                                                   deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, LALparams);
            break;

        case SEOBNRv2_opt:
             /* Waveform-specific sanity checks */
             if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                 XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
             if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                 XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
             if( !checkTidesZero(lambda1, lambda2) )
                 XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
             if( f_ref != 0.)
                 XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
             /* Call the waveform driver routine */
             SpinAlignedEOBversion = 200;
             ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                     deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, LALparams);
             break;

        case SEOBNRv3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
        REAL8 spin1[3], spin2[3];
            spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
            spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            PrecEOBversion = 3;
            ret = XLALSimIMRSpinEOBWaveform(hplus, hcross, /*&epoch,*/ phiRef,
                                            deltaT, m1, m2, f_min, distance, inclination, spin1, spin2, PrecEOBversion);
            break;

        case SEOBNRv4_opt:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 400;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                                                   deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, LALparams);
            break;

    case HGimri:
         /* Waveform-specific sanity checks */
         if( !checkTidesZero(lambda1, lambda2) )
         XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
         if( !checkCOSpinZero(S2x, S2y, S2z) )
         XLAL_ERROR(XLAL_EINVAL, "Non-zero CO spin given, but this approximant does not support this case.");
         /* Call the waveform driver */
         ret = XLALHGimriGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, S1z);
         break;

        case NR_hdf5:
            /* Waveform-specific sanity checks */
            numrel_data_path = XLALSimInspiralGetNumrelDataOLD(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralNRWaveformGetHplusHcross(hplus, hcross,
                    phiRef, inclination, deltaT, m1, m2, distance, f_min, f_ref, S1x, S1y, S1z,
                    S2x, S2y, S2z, numrel_data_path, NULL);
            XLALFree(numrel_data_path);
            break;


        default:
            XLALPrintError("TD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    if (polariz) {
      REAL8 tmpP,tmpC;
      REAL8 cp=cos(2.*polariz);
      REAL8 sp=sin(2.*polariz);
      for (UINT4 idx=0;idx<(*hplus)->data->length;idx++) {
    tmpP=(*hplus)->data->data[idx];
    tmpC=(*hcross)->data->data[idx];
    (*hplus)->data->data[idx] =cp*tmpP+sp*tmpC;
    (*hcross)->data->data[idx]=cp*tmpC-sp*tmpP;
      }
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

    return ret;
}

/**
 * @deprecated Use XLALSimInspiralChooseFDWaveform() instead
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the frequency domain.
 */
int XLALSimInspiralChooseFDWaveformOLD(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    const REAL8 m1,                         /**< mass of companion 1 (kg) */
    const REAL8 m2,                         /**< mass of companion 2 (kg) */
    const REAL8 S1x,                        /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                        /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                        /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                        /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                        /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                        /**< z-component of the dimensionless spin of object 2 */
    const REAL8 distance,                   /**< distance of source (m) */
    const REAL8 inclination,                /**< inclination of source (rad) */
    const REAL8 phiRef,                     /**< reference orbital phase (rad) */
    const REAL8 longAscNodes,               /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    const REAL8 eccentricity,               /**< eccentricity at reference epoch */
    const REAL8 UNUSED meanPerAno,          /**< mean anomaly of periastron */
    // frequency sampling parameters, no default value
    const REAL8 deltaF,                     /**< sampling interval (Hz) */
    const REAL8 f_min,                      /**< starting GW frequency (Hz) */
    const REAL8 f_max,                      /**< ending GW frequency (Hz) */
    REAL8 f_ref,                            /**< Reference frequency (Hz) */
    const REAL8 lambda1,                    /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    const REAL8 lambda2,                    /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    const REAL8 dQuadParam1,                /**< (quad-monop parameter of mass 1) / m1^5 -1 (dimensionless), give 0 for BHs */
    const REAL8 dQuadParam2,                /**< (quad-monop parameter of mass 2) / m2^5 -1 (dimensionless), give 0 for BHs */
    LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                         /**< twice post-Newtonian amplitude order */
    const int phaseO,                       /**< twice post-Newtonian order */
    const Approximant approximant           /**< post-Newtonian approximant to use for waveform production */
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
    REAL8 quadparam1 = 1.+dQuadParam1, quadparam2 = 1.+dQuadParam2;
    INT4 phiRefAtEnd;

    /* Support variables for precessing wfs*/
    REAL8 incl;
    REAL8 spin1x,spin1y,spin1z;
    REAL8 spin2x,spin2y,spin2z;

    /* Variables for IMRPhenomP and IMRPhenomPv2 */
    REAL8 chi1_l, chi2_l, chip, thetaJ, alpha0;

    /* General sanity checks that will abort
     *
     * If non-GR approximants are added, include them in
     * XLALSimInspiralApproximantAcceptTestGRParams()
     */
    if( nonGRparams && XLALSimInspiralApproximantAcceptTestGRParams(approximant) != LAL_SIM_INSPIRAL_TESTGR_PARAMS ) {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaF > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaF = %e requested...This corresponds to a very short TD signal (with padding). Consider a smaller value.\n", __func__, deltaF);
    if( deltaF < 1./4096. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaF = %e requested...This corresponds to a very long TD signal. Consider a larger value.\n", __func__, deltaF);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested...Signal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( S1x*S1x + S1y*S1y + S1z*S1z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if( S2x*S2x + S2y*S2y + S2z*S2z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested...Check for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
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
    pfac = 0.5 * (1. + cfac*cfac);

    switch (approximant)
    {
        /* inspiral-only models */
    case EccentricFD:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            /* Call the waveform driver routine */
            /* Note that for generic inclined eccentric waveforms
 *             * it is not possible to decompose hc(f) \propto I * hp(f)
 *                         * we call both polarizations independently
 *                                     */
            /*ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2,
 *             f_min, f_max, i, r, lambda1, lambda2, phaseO);*/
            ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2,
            f_min, f_max, inclination, distance,  (REAL8) XLALSimInspiralGetTestGRParam( nonGRparams, "inclination_azimuth"),
             eccentricity, phaseO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case TaylorF2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        XLAL_PRINT_DEPRECATION_WARNING("Calling TF2 via old interface, setting to default values tidal lambdas, quad-monopole pars, amplitude and phase order");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2(hptilde, phiRef, deltaF, m1, m2,
                    S1z, S2z, f_min, f_max, f_ref, distance,
                    NULL);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case TaylorF2NLTides:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            XLAL_PRINT_DEPRECATION_WARNING("Calling TF2 via old interface, setting to default values tidal lambdas, quad-monopole pars, amplitude and phase order");

            // FIXME : add checks for NL tidal parameters?

            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2NLTides(hptilde, phiRef, deltaF, m1, m2,
                    S1z, S2z, f_min, f_max, f_ref, distance,
                    NULL);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* non-spinning inspiral-merger-ringdown models */
        case IMRPhenomA:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomAGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    f_min, f_max, distance);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* spinning inspiral-only models */
        case SpinTaylorF2:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkCOSpinZero(S2x, S2y, S2z) ) // This is a single-spin model
                XLAL_ERROR(XLAL_EINVAL, "Non-zero CO spin given, but this approximant does not support this case.");
        spin1x=S1x;
        spin1y=S1y;
        spin1z=S1z;
            ROTATEY(inclination, spin1x, spin1y, spin1z);
            LNhatx = sin(inclination);
            LNhaty = 0.;
            LNhatz = cos(inclination);
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = 0; /* amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER */;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorF2(hptilde, hctilde, phiRef, deltaF,
                    m1, m2, spin1x, spin1y, spin1z, LNhatx, LNhaty, LNhatz,
                    f_min, f_max, f_ref, distance,
                    NULL, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        /* FIXME: Comment out this case, as I don't have its source code */
        //case TaylorR2F4:
        //    /* Waveform-specific sanity checks */
        //    if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
        //        XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
        //    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
        //        XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        //    /* Call the waveform driver routine */
        //    ret = XLALSimInspiralTaylorR2F4(hptilde, phiRef, deltaF, m1, m2,
        //            S1z, S2z, f_min, r, phaseO, amplitudeO);
        //    break;

        case TaylorF2RedSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpin(hptilde, phiRef, deltaF,
                    m1, m2, XLALSimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, distance, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case TaylorF2RedSpinTidal:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpinTidal(hptilde,phiRef,deltaF,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    lambda1, lambda2, f_min, f_max, distance, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* spinning inspiral-merger-ringdown models */
        case IMRPhenomB:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomBGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, distance);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case IMRPhenomC:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomCGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, distance, NULL);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case IMRPhenomD:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            /* Call the waveform driver routine */

            ret = XLALSimIMRPhenomDGenerateFD(hptilde, phiRef, f_ref, deltaF, m1, m2,
                  S1z, S2z, f_min, f_max, distance, NULL, NoNRT_V);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* Produce both polarizations */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case EOBNRv2_ROM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

            ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde,
        phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, 0);
            break;

        case EOBNRv2HM_ROM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

            ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde,
        phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, 1);
            break;

        case SEOBNRv1_ROM_EffectiveSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if (!checkAlignedSpinsEqual(S1z, S2z)) {
                    XLALPrintError("XLAL Error - %s: SEOBNRv1ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__,S1z,S2z);
                    XLAL_ERROR(XLAL_EINVAL);
            }

            ret = XLALSimIMRSEOBNRv1ROMEffectiveSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
            break;

        case SEOBNRv1_ROM_DoubleSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

            ret = XLALSimIMRSEOBNRv1ROMDoubleSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z);
            break;

        case SEOBNRv2_ROM_EffectiveSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            if (!checkAlignedSpinsEqual(S1z, S2z)) {
                    XLALPrintError("XLAL Error - %s: SEOBNRv2ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__,S1z,S2z);
                    XLAL_ERROR(XLAL_EINVAL);
            }

            ret = XLALSimIMRSEOBNRv2ROMEffectiveSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
            break;

        case SEOBNRv2_ROM_DoubleSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

            ret = XLALSimIMRSEOBNRv2ROMDoubleSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z);
            break;

        case SEOBNRv2_ROM_DoubleSpin_HI:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefaultOLD(waveFlags) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.");
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

            ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHI(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1);
            break;


        case IMRPhenomP:
            /* Waveform-specific sanity checks */
            XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
            if( !XLALSimInspiralModesChoiceIsDefault(          /* Default is (2,2) or l=2 modes. */
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            LNhatx = sin(incl);
            LNhaty = 0.;
            LNhatz = cos(incl);
            /* Tranform to model parameters */
            if(f_ref==0.0)
                f_ref = f_min; /* Default reference frequency is minimum frequency */
            XLALSimIMRPhenomPCalculateModelParametersOld(
                &chi1_l, &chi2_l, &chip, &thetaJ, &alpha0,
                m1, m2, f_ref,
                LNhatx, LNhaty, LNhatz,
        spin1x, spin1y, spin1z,
        spin2x, spin2y, spin2z, IMRPhenomPv1_V);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomP(hptilde, hctilde,
              chi1_l, chi2_l, chip, thetaJ,
              m1, m2, distance, alpha0, phiRef, deltaF, f_min, f_max, f_ref, IMRPhenomPv1_V, NoNRT_V, NULL);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case IMRPhenomPv2:
            /* Waveform-specific sanity checks */
      XLALSimInspiralInitialConditionsPrecessingApproxs(&incl,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,inclination,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,f_ref,phiRef,XLALSimInspiralGetFrameAxis(waveFlags));
            if( !XLALSimInspiralModesChoiceIsDefault(          /* Default is (2,2) or l=2 modes. */
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
            if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
            LNhatx = sin(incl);
            LNhaty = 0.;
            LNhatz = cos(incl);
            /* Tranform to model parameters */
            if(f_ref==0.0)
                f_ref = f_min; /* Default reference frequency is minimum frequency */
            XLALSimIMRPhenomPCalculateModelParametersOld(
                &chi1_l, &chi2_l, &chip, &thetaJ, &alpha0,
                m1, m2, f_ref,
                LNhatx, LNhaty, LNhatz,
                spin1x, spin1y, spin1z,
                spin2x, spin2y, spin2z, IMRPhenomPv2_V);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomP(hptilde, hctilde,
              chi1_l, chi2_l, chip, thetaJ,
              m1, m2, distance, alpha0, phiRef, deltaF, f_min, f_max, f_ref, IMRPhenomPv2_V, NoNRT_V, NULL);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case SpinTaylorT4Fourier:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        spin1x=S1x;
        spin1y=S1y;
        spin1z=S1z;
        spin2x=S2x;
        spin2y=S2y;
        spin2z=S2z;
        ROTATEY(inclination,spin1x,spin1y,spin1z);
        ROTATEY(inclination,spin2x,spin2y,spin2z);
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
            fStart = 0.9*f_min;
            phiRefAtEnd = 0;
            // if f_ref = 0, set it to f_min, and tell the driver routine that we came from there
            if(f_ref == 0)
            {
              f_ref = f_min;
              phiRefAtEnd = 1;
            }
            // default quadparams are for black holes. Replace by ~2-12 for neutron stars
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT4Fourier(hptilde, hctilde,
              f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2, NULL, phaseO, amplitudeO, phiRefAtEnd);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case SpinTaylorT5Fourier:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        spin1x=S1x;
        spin1y=S1y;
        spin1z=S1z;
        spin2x=S2x;
        spin2y=S2y;
        spin2z=S2z;
        ROTATEY(inclination,spin1x,spin1y,spin1z);
        ROTATEY(inclination,spin2x,spin2y,spin2z);
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
            fStart = 0.9*f_min;
            phiRefAtEnd = 0;
            // if f_ref = 0, set it to f_min, and tell the driver routine that we came from there
            if(f_ref == 0)
            {
              f_ref = f_min;
              phiRefAtEnd = 1;
            }
            // default quadparams are for black holes. Replace by ~2-12 for neutron stars
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT5Fourier(hptilde, hctilde,
                             f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2, NULL, phaseO, amplitudeO, phiRefAtEnd);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;


        default:
            XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 polariz=longAscNodes;
    if (polariz) {
      COMPLEX16 tmpP,tmpC;
      for (UINT4 idx=0;idx<(*hptilde)->data->length;idx++) {
    tmpP=(*hptilde)->data->data[idx];
    tmpC=(*hctilde)->data->data[idx];
    (*hptilde)->data->data[idx] =cos(2.*polariz)*tmpP+sin(2.*polariz)*tmpC;
    (*hctilde)->data->data[idx]=cos(2.*polariz)*tmpC-sin(2.*polariz)*tmpP;
      }
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

    return ret;
}

/**
 * if you do NOT provide a quadparam[1,2] term and you DO provide
 * lamdba[1,2] then we calculate quad-mono term using universal relations
 * quadparam[1,2]_UR: Quadrupole-Monopole parameter computed using
 * universal relations (UR)
 */

int XLALSimInspiralSetQuadMonParamsFromLambdas(
       LALDict *LALparams /**< LAL dictionary containing accessory parameters */
       )
{
    REAL8 quadparam1 = XLALSimInspiralWaveformParamsLookupdQuadMon1(LALparams);
    REAL8 quadparam2 = XLALSimInspiralWaveformParamsLookupdQuadMon2(LALparams);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALparams);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALparams);

    if ((lambda1 > 0) && (quadparam1 == 0)) {
        REAL8 quadparam1_UR = XLALSimInspiralEOSQfromLambda(lambda1);
        XLALSimInspiralWaveformParamsInsertdQuadMon1(LALparams, quadparam1_UR - 1.);
    }

    if ((lambda2 > 0) && (quadparam2 == 0)) {
        REAL8 quadparam2_UR = XLALSimInspiralEOSQfromLambda(lambda2);
        XLALSimInspiralWaveformParamsInsertdQuadMon2(LALparams, quadparam2_UR - 1.);
    }
    return XLAL_SUCCESS;
}
/** @} */
