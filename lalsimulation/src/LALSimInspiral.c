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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
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

#include "LALSimInspiralPNCoefficients.c"
#include "check_series_macros.h"
#include "check_waveform_macros.h"

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

/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) - vy*sin(angle);\
	tmp2 = vx*sin(angle) + vy*cos(angle);\
	vx = tmp1;\
	vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) + vz*sin(angle);\
	tmp2 = - vx*sin(angle) + vz*cos(angle);\
	vx = tmp1;\
	vz = tmp2


#define INITIALIZE_NAME(a) [a] = #a
/* TODO: UPDATE ME WHENEVER A NEW APPROXIMANT IS ADDED */
static const char *lalSimulationApproximantNames[] = {
    INITIALIZE_NAME(TaylorT1),
    INITIALIZE_NAME(TaylorT2),
    INITIALIZE_NAME(TaylorT3),
    INITIALIZE_NAME(TaylorF1),
    INITIALIZE_NAME(TaylorF2),
    INITIALIZE_NAME(TaylorR2F4),
    INITIALIZE_NAME(TaylorF2RedSpin),
    INITIALIZE_NAME(TaylorF2RedSpinTidal),
    INITIALIZE_NAME(PadeT1),
    INITIALIZE_NAME(PadeF1),
    INITIALIZE_NAME(EOB),
    INITIALIZE_NAME(BCV),
    INITIALIZE_NAME(BCVSpin),
    INITIALIZE_NAME(SpinTaylorT1),
    INITIALIZE_NAME(SpinTaylorT2),
    INITIALIZE_NAME(SpinTaylorT3),
    INITIALIZE_NAME(SpinTaylorT4),
    INITIALIZE_NAME(SpinTaylorT5),
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
    INITIALIZE_NAME(SEOBNRv1),
    INITIALIZE_NAME(SEOBNRv2),
    INITIALIZE_NAME(SEOBNRv2_opt),
    INITIALIZE_NAME(SEOBNRv3),
    INITIALIZE_NAME(SEOBNRv4),
    INITIALIZE_NAME(SEOBNRv4_opt),
    INITIALIZE_NAME(SEOBNRv1_ROM_EffectiveSpin),
    INITIALIZE_NAME(SEOBNRv1_ROM_DoubleSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_EffectiveSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_DoubleSpin),
    INITIALIZE_NAME(SEOBNRv2_ROM_DoubleSpin_HI),
    INITIALIZE_NAME(HGimri),
    INITIALIZE_NAME(IMRPhenomA),
    INITIALIZE_NAME(IMRPhenomB),
    INITIALIZE_NAME(IMRPhenomFA),
    INITIALIZE_NAME(IMRPhenomFB),
    INITIALIZE_NAME(IMRPhenomC),
    INITIALIZE_NAME(IMRPhenomD),
    INITIALIZE_NAME(IMRPhenomP),
    INITIALIZE_NAME(IMRPhenomPv2),
    INITIALIZE_NAME(IMRPhenomFC),
    INITIALIZE_NAME(TaylorEt),
    INITIALIZE_NAME(TaylorT4),
    INITIALIZE_NAME(EccentricTD),
    INITIALIZE_NAME(TaylorN),
    INITIALIZE_NAME(SpinTaylorT4Fourier),
    INITIALIZE_NAME(SpinTaylorT2Fourier),
    INITIALIZE_NAME(SpinDominatedWf),
    INITIALIZE_NAME(NR_hdf5),
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


/*
 * certain approximants adopt the convention that f_ref=0 refers to the start
 * of the waveform while other approximants adopt the convention that f_ref=0
 * refers to the end of the waveform. in the former case, this routine will
 * return the explicit value of f_ref, which will be f_min.
 */
static double fixReferenceFrequency(double f_ref, double f_min, Approximant approximant)
{
    if (f_ref == 0)
        switch (approximant) {
        case SpinTaylorT1:
        case SpinTaylorT2:
        case SpinTaylorT3:
        case SpinTaylorT4:
        case SpinTaylorT2Fourier:
        case SpinTaylorT4Fourier:
        case SpinTaylorF2:
        case IMRPhenomP:
        case IMRPhenomPv2:
            return f_min;
        default:
            break;
        }
    return f_ref;
}

/*
 * some helper routines for XLALSimInspiralTD
 */
static int XLALSimInspiralTDFromTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 z, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
static int XLALSimInspiralTDFromFD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 z, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);


/**
 * @addtogroup LALSimInspiral_c
 * @brief General routines for generating binary inspiral waveforms.
 *
 * @{
 */

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
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 LNhatx, LNhaty, LNhatz, E1x, E1y, E1z;
    char* numrel_data_path;
    //REAL8 tmp1, tmp2;
    int ret;
    /* N.B. the quadrupole of a spinning compact body labeled by A is
     * Q_A = - quadparam_A chi_A^2 m_A^3 (see gr-qc/9709032)
     * where quadparam = 1 for BH ~= 4-8 for NS.
     * This affects the quadrupole-monopole interaction.
     * For now, hardcode quadparam1,2 = 1.
     * Will later add ability to set via LALSimInspiralTestGRParam
     */
    REAL8 v0 = 1., quadparam1 = 1., quadparam2 = 1.;

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
    REAL8 iTmp;

    /* SEOBNR flag for model version. 1 for SEOBNRv1, 2 for SEOBNRv2 */
    UINT4 SpinAlignedEOBversion;
    REAL8 spin1[3], spin2[3];
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
    f_ref = fixReferenceFrequency(f_ref, f_min, approximant);

    switch (approximant)
    {
        /* non-spinning inspiral-only models */
        case TaylorEt:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, r, i, amplitudeO, phaseO);
            break;

        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT2PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

	case EccentricTD:
	    /* Waveform-specific sanity checks */
	    if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
    		ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
	    if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
    		ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
	    if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
    	    	ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
	    if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
    		ABORT_NONZERO_SPINS(waveFlags);
	    if( !checkTidesZero(lambda1, lambda2) )
		ABORT_NONZERO_TIDES(waveFlags);
	    /* Call the waveform driver routine */
	    ret = XLALSimInspiralEccentricTDPNGenerator(hplus, hcross, phiRef,
		    deltaT, m1, m2, f_min, f_ref, r, i, (REAL8) XLALSimInspiralGetTestGRParam( nonGRparams, "e_min"),
		    amplitudeO, phaseO);
	    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
	    break;

        /* non-spinning inspiral-merger-ringdown models */
        case IMRPhenomA:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomAGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, 0., r, i);
            break;

        case EOBNRv2HM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // FIXME: need to create a function to take in different modes or produce an error if all modes not given
            ret = XLALSimIMREOBNRv2AllModes(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, r, i);
            break;

        case EOBNRv2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMREOBNRv2DominantMode(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, r, i);
            break;

        /* spinning inspiral-only models */
        case SpinTaylorT2:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
	    spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
	    spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
	    iTmp=i;
	    XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = - sin(i);
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT2(hplus, hcross, phiRef, v0, deltaT,
                    m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO);
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
	    spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
	    spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
	    iTmp=i;
	    XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = - sin(i);
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT4(hplus, hcross, phiRef, v0, deltaT,
                    m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO);
            break;

	 case SpinTaylorT1:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
	    spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
	    spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
	    iTmp=i;
	    XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = - sin(i);
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT1(hplus, hcross, phiRef, v0, deltaT,
                    m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
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
		spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
	        spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
	        iTmp=i;
	        XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
                LNhatx = sin(i);
                LNhaty = 0.;
                LNhatz = cos(i);
                /* Call the waveform driver routine */
                ret = XLALSimInspiralSpinDominatedWaveformInterfaceTD(hplus, hcross, deltaT, m1, m2, f_min, f_ref, r, S1x, S1y, S1z, LNhatx, LNhaty, LNhatz, phaseO, amplitudeO, phiRef);
                break;

        /* spin aligned inspiral-merger-ringdown models */
        case IMRPhenomB:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomBGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, 0., r, i);
            break;

        case PhenSpinTaylor:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimSpinInspiralGenerator(hplus, hcross, phiRef,
					       deltaT, m1, m2, f_min, f_ref, r, i, S1x, S1y, S1z, S2x, S2y, S2z,
					       phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, waveFlags, nonGRparams);
	    break;

        case IMRPhenomC:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomCGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, 0., r, i, nonGRparams);
            break;

	case IMRPhenomD:
	    if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
		    ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
	    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
		    ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
	    if( !checkTidesZero(lambda1, lambda2) )
		    ABORT_NONZERO_TIDES(waveFlags);
	    // generate TD waveforms with zero inclincation so that amplitude can be
	    // calculated from hplus and hcross, apply inclination-dependent factors
	    // in loop below
	    ret = XLALSimInspiralTDFromFD(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z,
			    S2x, S2y, S2z, f_min, f_ref, r, 0, 0, lambda1, lambda2,
			    waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
	    REAL8 maxamp=0;
	    REAL8TimeSeries *hp = *hplus;
	    REAL8TimeSeries *hc = *hcross;
	    INT4 maxind=hp->data->length - 1;
	    INT4 loopi;
	    const REAL8 cfac=cos(i);
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
	    ret = XLALSimInspiralTDFromFD(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z,
			    S2x, S2y, S2z, f_min, f_ref, r, 0, i, lambda1, lambda2,
			    waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
	    break;

        case PhenSpinTaylorRD:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at the start.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenSpinInspiralRDGenerator(hplus, hcross, phiRef,
							deltaT, m1, m2, f_min, f_ref, r, i, S1x, S1y, S1z, S2x, S2y, S2z,
							phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, waveFlags, nonGRparams);
            break;

        case SEOBNRv1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 1;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, f_min, r, i, S1z, S2z, SpinAlignedEOBversion);
            break;

        case SEOBNRv2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 2;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, f_min, r, i, S1z, S2z, SpinAlignedEOBversion);
            break;

        case SEOBNRv4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 4;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                                                   deltaT, m1, m2, f_min, r, i, S1z, S2z, SpinAlignedEOBversion);
            break;

        case SEOBNRv2_opt:
             /* Waveform-specific sanity checks */
             if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                 ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
             if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                 ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
             if( !checkTidesZero(lambda1, lambda2) )
                 ABORT_NONZERO_TIDES(waveFlags);
             if( f_ref != 0.)
                 XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
             /* Call the waveform driver routine */
             SpinAlignedEOBversion = 200;
             ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                     deltaT, m1, m2, f_min, r, i, S1z, S2z, SpinAlignedEOBversion);
             break;

        case SEOBNRv3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
            spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            iTmp=i;
           //XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
            spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            ret = XLALSimIMRSpinEOBWaveform(hplus, hcross, /*&epoch,*/ phiRef,
                    deltaT, m1, m2, f_min, r, i, spin1, spin2);
            break;

        case SEOBNRv4_opt:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            SpinAlignedEOBversion = 400;
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef,
                                                   deltaT, m1, m2, f_min, r, i, S1z, S2z, SpinAlignedEOBversion);
            break;

	case HGimri:
	     /* Waveform-specific sanity checks */
	     if( !checkTidesZero(lambda1, lambda2) )
		 ABORT_NONZERO_TIDES(waveFlags);
	     if( !checkCOSpinZero(S2x, S2y, S2z) )
		 ABORT_NONZERO_SPIN2(waveFlags);
	     /* Call the waveform driver */
	     ret = XLALHGimriGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, r, i, S1z);
	     break;

        case NR_hdf5:
            /* Waveform-specific sanity checks */

            /* Call the waveform driver routine */
            numrel_data_path = XLALSimInspiralGetNumrelData(waveFlags);
            ret = XLALSimInspiralNRWaveformGetHplusHcross(hplus, hcross,
                    phiRef, i, deltaT, m1, m2, r, f_min, f_ref, S1x, S1y, S1z,
                    S2x, S2y, S2z, numrel_data_path);
            XLALFree(numrel_data_path);
            break;


        default:
            XLALPrintError("TD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

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
    REAL8 phiRef,                           /**< reference orbital phase (rad) */
    REAL8 deltaF,                           /**< sampling interval (Hz) */
    REAL8 m1,                               /**< mass of companion 1 (kg) */
    REAL8 m2,                               /**< mass of companion 2 (kg) */
    REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                            /**< starting GW frequency (Hz) */
    REAL8 f_max,                            /**< ending GW frequency (Hz) */
    REAL8 f_ref,                            /**< Reference frequency (Hz) */
    REAL8 r,                                /**< distance of source (m) */
    REAL8 i,                                /**< inclination of source (rad) */
    REAL8 lambda1,                          /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                          /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                         /**< twice post-Newtonian amplitude order */
    int phaseO,                             /**< twice post-Newtonian order */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
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
    REAL8 quadparam1 = 1., quadparam2 = 1.; /* FIXME: This cannot yet be set in the interface */
    INT4 phiRefAtEnd;

    /* Support variables for precessing wfs*/
    REAL8 iTmp;
    REAL8 spin1[3],spin2[3];

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
    f_ref = fixReferenceFrequency(f_ref, f_min, approximant);

    /* The non-precessing waveforms return h(f) for optimal orientation
     * (i=0, Fp=1, Fc=0; Lhat pointed toward the observer)
     * To get generic polarizations we multiply by inclination dependence
     * and note hc(f) \propto -I * hp(f)
     * Non-precessing waveforms multiply hp by pfac, hc by -I*cfac
     */
    cfac = cos(i);
    pfac = 0.5 * (1. + cfac*cfac);

    switch (approximant)
    {
        /* inspiral-only models */
	case EccentricFD:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            /* Call the waveform driver routine */
            /* Note that for generic inclined eccentric waveforms
 *             * it is not possible to decompose hc(f) \propto I * hp(f)
 *                         * we call both polarizations independently
 *                                     */
            /*ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2,
 *             f_min, f_max, i, r, lambda1, lambda2, phaseO);*/
            ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2,
            f_min, f_max, i, r,  (REAL8) XLALSimInspiralGetTestGRParam( nonGRparams, "inclination_azimuth"),
             (REAL8) XLALSimInspiralGetTestGRParam( nonGRparams, "e_min"),
              phaseO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case TaylorF2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);

            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2(hptilde, phiRef, deltaF, m1, m2,
                    S1z, S2z, f_min, f_max, f_ref, r,
                    quadparam1, quadparam2, lambda1, lambda2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO, nonGRparams);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomAGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    f_min, f_max, r);
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
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkCOSpinZero(S2x, S2y, S2z) ) // This is a single-spin model
                ABORT_NONZERO_SPIN2(waveFlags);
            ROTATEY(i, S1x, S1y, S1z);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            /* Maximum PN amplitude order for precessing waveforms is
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = 0; /* amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER */;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorF2(hptilde, hctilde, phiRef, deltaF,
                    m1, m2, S1x, S1y, S1z, LNhatx, LNhaty, LNhatz,
                    f_min, f_max, f_ref, r,
                    nonGRparams, XLALSimInspiralGetSpinOrder(waveFlags), phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        /* FIXME: Comment out this case, as I don't have its source code */
        //case TaylorR2F4:
        //    /* Waveform-specific sanity checks */
        //    if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
        //        ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
        //    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
        //        ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
        //    /* Call the waveform driver routine */
        //    ret = XLALSimInspiralTaylorR2F4(hptilde, phiRef, deltaF, m1, m2,
        //            S1z, S2z, f_min, r, phaseO, amplitudeO);
        //    break;

        case TaylorF2RedSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpin(hptilde, phiRef, deltaF,
                    m1, m2, XLALSimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r, phaseO, amplitudeO);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpinTidal(hptilde,phiRef,deltaF,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    lambda1, lambda2, f_min, f_max, r, phaseO, amplitudeO);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomBGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomCGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r, nonGRparams);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomDGenerateFD(hptilde, phiRef, f_ref, deltaF, m1, m2,
                  S1z, S2z, f_min, f_max, r, nonGRparams);
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
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);

            ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde,
		phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, 0);
            break;

        case EOBNRv2HM_ROM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);

            ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde,
		phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, 1);
            break;

        case SEOBNRv1_ROM_EffectiveSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if (!checkAlignedSpinsEqual(S1z, S2z)) {
                    XLALPrintError("XLAL Error - %s: SEOBNRv1ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__,S1z,S2z);
                    XLAL_ERROR(XLAL_EINVAL);
            }

            ret = XLALSimIMRSEOBNRv1ROMEffectiveSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
            break;

        case SEOBNRv1_ROM_DoubleSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);

            ret = XLALSimIMRSEOBNRv1ROMDoubleSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, S1z, S2z);
            break;

        case SEOBNRv2_ROM_EffectiveSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if (!checkAlignedSpinsEqual(S1z, S2z)) {
                    XLALPrintError("XLAL Error - %s: SEOBNRv2ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__,S1z,S2z);
                    XLAL_ERROR(XLAL_EINVAL);
            }

            ret = XLALSimIMRSEOBNRv2ROMEffectiveSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
            break;

        case SEOBNRv2_ROM_DoubleSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);

            ret = XLALSimIMRSEOBNRv2ROMDoubleSpin(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, S1z, S2z);
            break;

        case SEOBNRv2_ROM_DoubleSpin_HI:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);

            ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHI(hptilde, hctilde,
                    phiRef, deltaF, f_min, f_max, f_ref, r, i, m1, m2, S1z, S2z, -1);
            break;


        case IMRPhenomP:
            /* Waveform-specific sanity checks */
            spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
            spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
            iTmp=i;
            XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            if( !XLALSimInspiralModesChoiceIsDefault(          /* Default is (2,2) or l=2 modes. */
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            /* Tranform to model parameters */
            if(f_ref==0.0)
                f_ref = f_min; /* Default reference frequency is minimum frequency */
            XLALSimIMRPhenomPCalculateModelParameters(
                &chi1_l, &chi2_l, &chip, &thetaJ, &alpha0,
                m1, m2, f_ref,
                LNhatx, LNhaty, LNhatz,
                S1x, S1y, S1z,
                S2x, S2y, S2z, IMRPhenomPv1_V);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomP(hptilde, hctilde,
              chi1_l, chi2_l, chip, thetaJ,
              m1, m2, r, alpha0, phiRef, deltaF, f_min, f_max, f_ref, IMRPhenomPv1_V, nonGRparams);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case IMRPhenomPv2:
            /* Waveform-specific sanity checks */
            spin1[0]=S1x; spin1[1]=S1y; spin1[2]=S1z;
            spin2[0]=S2x; spin2[1]=S2y; spin2[2]=S2z;
            iTmp=i;
            XLALSimInspiralInitialConditionsPrecessingApproxs(&i,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,iTmp,spin1[0],spin1[1],spin1[2],spin2[0],spin2[1],spin2[2],m1,m2,f_ref,XLALSimInspiralGetFrameAxis(waveFlags));
            if( !XLALSimInspiralModesChoiceIsDefault(          /* Default is (2,2) or l=2 modes. */
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            /* Tranform to model parameters */
            if(f_ref==0.0)
                f_ref = f_min; /* Default reference frequency is minimum frequency */
            XLALSimIMRPhenomPCalculateModelParameters(
                &chi1_l, &chi2_l, &chip, &thetaJ, &alpha0,
                m1, m2, f_ref,
                LNhatx, LNhaty, LNhatz,
                S1x, S1y, S1z,
                S2x, S2y, S2z, IMRPhenomPv2_V);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomP(hptilde, hctilde,
              chi1_l, chi2_l, chip, thetaJ,
              m1, m2, r, alpha0, phiRef, deltaF, f_min, f_max, f_ref, IMRPhenomPv2_V, nonGRparams);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case SpinTaylorT4Fourier:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
	    ROTATEY(i,S1x,S1y,S1z);
	    ROTATEY(i,S2x,S2y,S2z);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = -sin(i);
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
            quadparam1 = 1.;
            quadparam2 = 1.;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT4Fourier(hptilde, hctilde,
              f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2, XLALSimInspiralGetSpinOrder(waveFlags), XLALSimInspiralGetTidalOrder(waveFlags), phaseO, amplitudeO, phiRefAtEnd);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;

        case SpinTaylorT2Fourier:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
	    ROTATEY(i,S1x,S1y,S1z);
	    ROTATEY(i,S2x,S2y,S2z);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = -sin(i);
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
            quadparam1 = 1.;
            quadparam2 = 1.;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT2Fourier(hptilde, hctilde,
              f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2, XLALSimInspiralGetSpinOrder(waveFlags), XLALSimInspiralGetTidalOrder(waveFlags), phaseO, amplitudeO, phiRefAtEnd);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            break;


        default:
            XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

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
 * This routine has one additional parameter relative to XLALSimInspiralChooseTDWaveform.
 * The additional parameter is the redshift, z, of the waveform.  This should be set to
 * zero for sources in the nearby universe (that are nearly at rest relative to the
 * earth).  For sources at cosmological distances, the mass parameters m1 and m2 should
 * be interpreted as the physical (source frame) masses of the bodies and the distance
 * parameter r is the comoving (transverse) distance.  If the calling routine has already
 * applied cosmological "corrections" to m1 and m2 and regards r as a luminosity distance
 * then the redshift factor should again be set to zero.
 *
 * @note The parameters passed must be in SI units.
 */
SphHarmTimeSeries* XLALSimInspiralTDModesFromPolarizations(
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 z,                                    /**< redshift of source frame relative to observer frame */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8TimeSeries *hplus, *hcross;
    COMPLEX16TimeSeries *h22,*h2m2;
    SphHarmTimeSeries *hlm;
    UINT4 j;
    int retval;
    float fac = XLALSpinWeightedSphericalHarmonic(0., 0., -2, 2,2);

    /* Generate waveform via on-axis emission. Assumes only (2,2) and (2,-2) emission */
    retval = XLALSimInspiralTD(&hplus, &hcross, 0.0, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, z, 0.0, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
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
static int XLALSimInspiralTDFromTD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 z,                                    /**< redshift of source frame relative to observer frame */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
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
    f_ref = fixReferenceFrequency(f_ref, f_min, approximant);

    /* apply redshift correction to dimensionful source-frame quantities */
    if (z != 0.0) {
        m1 *= (1.0 + z);
        m2 *= (1.0 + z);
        r *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
    }
    /* set redshift to zero so we don't accidentally apply it again later */
    z = 0.0;

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
    retval = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, fstart, f_ref, r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
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
static int XLALSimInspiralTDFromFD(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 z,                                    /**< redshift of source frame relative to observer frame */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
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
    double tchirp, tmerge, textra;
    double fisco, fstart;
    double s;
    int retval;

    if (!XLALSimInspiralImplementedFDApproximants(approximant))
        XLAL_ERROR(XLAL_EINVAL, "Invalid approximant: not a FD approximant");

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    f_ref = fixReferenceFrequency(f_ref, f_min, approximant);

    /* apply redshift correction to dimensionful source-frame quantities */
    if (z != 0.0) {
        m1 *= (1.0 + z);
        m2 *= (1.0 + z);
        r *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
    }
    /* set redshift to zero so we don't accidentally apply it again later */
    z = 0.0;

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
    retval = XLALSimInspiralFD(&hptilde, &hctilde, phiRef, 0.0, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, 0.5/deltaT, f_ref, r, 0.0, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
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
 *
 * This routine has one additional parameter relative to XLALSimInspiralChooseTDWaveform.
 * The additional parameter is the redshift, z, of the waveform.  This should be set to
 * zero for sources in the nearby universe (that are nearly at rest relative to the
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
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 z,                                    /**< redshift of source frame relative to observer frame */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    /* call the appropriate helper routine */
    if (XLALSimInspiralImplementedTDApproximants(approximant)) {
        if (XLALSimInspiralTDFromTD(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, z, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant) < 0)
            XLAL_ERROR(XLAL_EFUNC);
    } else if (XLALSimInspiralImplementedFDApproximants(approximant)) {
        if (XLALSimInspiralTDFromFD(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, z, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant) < 0)
            XLAL_ERROR(XLAL_EFUNC);
    } else
        XLAL_ERROR(XLAL_EINVAL, "Invalid approximant");
    return 0;
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
 * a deltaF that is small enough to represent the Fourier transform of a time-domain waveform.
 * If deltaF is specified but f_max / deltaF is not a power of 2, and the waveform approximant
 * is a time-domain approximant, then f_max is increased so that f_max / deltaF is the next
 * power of 2.  (If the user wishes to discard the extra high frequency content, this must
 * be done separately.)
 *
 * This routine has one additional parameter relative to XLALSimInspiralChooseFDWaveform.
 * The additional parameter is the redshift, z, of the waveform.  This should be set to
 * zero for sources in the nearby universe (that are nearly at rest relative to the
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
    COMPLEX16FrequencySeries **hptilde,         /**< +-polarization waveform */
    COMPLEX16FrequencySeries **hctilde,         /**< x-polarization waveform */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaF,                               /**< frequency interval (Hz), or 0 to compute necessary interval */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_max,                                /**< ending GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 z,                                    /**< redshift of source frame relative to observer frame */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double chirplen, deltaT;
    int chirplen_exp;
    int retval;

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    f_ref = fixReferenceFrequency(f_ref, f_min, approximant);

    /* apply redshift correction to dimensionful source-frame quantities */
    if (z != 0.0) {
        m1 *= (1.0 + z);
        m2 *= (1.0 + z);
        r *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
    }
    /* set redshift to zero so we don't accidentally apply it again later */
    z = 0.0;

    /* FIXME: assume that f_max is the Nyquist frequency, and use it
     * to compute the requested deltaT */
    deltaT = 0.5 / f_max;

    if (XLALSimInspiralImplementedFDApproximants(approximant)) {

        /* generate a FD waveform and condition it by applying tapers at
         * frequencies between a frequency below the requested f_min and
         * f_min; also wind the waveform in phase in case it would wrap-
         * around at the merger time */

        double tchirp, tmerge, textra, tshift;
        double fstart, fisco;
        double s;
        size_t k, k0, k1;

        /* if the requested low frequency is below the lowest Kerr ISCO
         * frequency then change it to that frequency */
        fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
        if (f_min > fisco)
            f_min = fisco;

        /* upper bound on the chirp time starting at f_min */
        tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);

        /* upper bound on the final plunge, merger, and ringdown time */
        switch (approximant) {
        case TaylorF2:
        case SpinTaylorF2:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
        case SpinTaylorT4Fourier:
            /* inspiral-only models: no merger time */
            tmerge = 0.0;
            break;
        default:
            /* IMR model: estimate plunge and merger time */
            /* sometimes these waveforms have phases that
             * cause them to wrap-around an amount equal to
             * the merger-ringodwn time, so we will undo
             * that here */
            s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);
            tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);
            break;
        }

        /* new lower frequency to start the waveform: add some extra early
         * part over which tapers may be applied, the extra amount being
         * a fixed fraction of the chirp time; add some additional padding
         * equal to a few extra cycles at the low frequency as well for
         * safety and for other routines to use */
        textra = extra_cycles / f_min;
        fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);

        /* revise (over-)estimate of chirp from new start frequency */
        tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, S1z, S2z);

        /* need a long enough segment to hold a whole chirp with some padding */
        /* length of the chirp in samples */
        chirplen = round((tchirp + tmerge + 2.0 * textra) / deltaT);
        /* make chirplen next power of two */
        frexp(chirplen, &chirplen_exp);
        chirplen = ldexp(1.0, chirplen_exp);
        /* frequency resolution */
        if (deltaF == 0.0)
            deltaF = 1.0 / (chirplen * deltaT);
        else if (deltaF > 1.0 / (chirplen * deltaT))
            XLAL_PRINT_WARNING("Specified frequency interval of %g Hz is too large for a chirp of duration %g s", deltaF, chirplen * deltaT);

        /* generate the waveform in the frequency domain starting at fstart */
        retval = XLALSimInspiralChooseFDWaveform(hptilde, hctilde, phiRef, deltaF, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, fstart, f_max, f_ref, r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
        if (retval < 0)
            XLAL_ERROR(XLAL_EFUNC);

        /* taper frequencies between fstart and f_min */
        k0 = round(fstart / (*hptilde)->deltaF);
        k1 = round(f_min / (*hptilde)->deltaF);
        /* make sure it is zero below fstart */
        for (k = 0; k < k0; ++k) {
            (*hptilde)->data->data[k] = 0.0;
            (*hctilde)->data->data[k] = 0.0;
        }
        /* taper between fstart and f_min */
        for ( ; k < k1; ++k) {
            double w = 0.5 - 0.5 * cos(M_PI * (k - k0) / (double)(k1 - k0));
            (*hptilde)->data->data[k] *= w;
            (*hctilde)->data->data[k] *= w;
        }
        /* make sure Nyquist frequency is zero */
        (*hptilde)->data->data[(*hptilde)->data->length - 1] = 0.0;
        (*hctilde)->data->data[(*hctilde)->data->length - 1] = 0.0;

        /* we want to make sure that this waveform will give something
         * sensible if it is later transformed into the time domain:
         * to avoid the end of the waveform wrapping around to the beginning,
         * we shift waveform backwards in time and compensate for this
         * shift by adjusting the epoch */
        tshift = round(tmerge / deltaT) * deltaT; /* integer number of time samples */
        for (k = 0; k < (*hptilde)->data->length; ++k) {
            double complex phasefac = cexp(2.0 * M_PI * I * k * deltaF * tshift);
            (*hptilde)->data->data[k] *= phasefac;
            (*hctilde)->data->data[k] *= phasefac;
        }
        XLALGPSAdd(&(*hptilde)->epoch, tshift);
        XLALGPSAdd(&(*hctilde)->epoch, tshift);

    } else if (XLALSimInspiralImplementedTDApproximants(approximant)) {

        /* generate a conditioned waveform in the time domain and Fourier-transform it */

        REAL8TimeSeries *hplus = NULL;
        REAL8TimeSeries *hcross = NULL;
        REAL8FFTPlan *plan;

        /* generate conditioned waveform in time domain */
        /* note: redshift factor has already been applied */
        retval = XLALSimInspiralTD(&hplus, &hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, 0.0, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
        if (retval < 0)
            XLAL_ERROR(XLAL_EFUNC);

        /* determine chirp length and round up to next power of two */
        chirplen = hplus->data->length;
        frexp(chirplen, &chirplen_exp);
        chirplen = ldexp(1.0, chirplen_exp);
        /* frequency resolution */
        if (deltaF == 0.0)
            deltaF = 1.0 / (chirplen * hplus->deltaT);
        else { /* recompute chirplen based on deltaF and f_max */
            size_t n;
            if (deltaF > 1.0 / (chirplen * deltaT))
                XLAL_PRINT_WARNING("Specified frequency interval of %g Hz is too large for a chirp of duration %g s", deltaF, chirplen * deltaT);
            n = chirplen = round(2.0 * f_max / deltaF);
            if ((n & (n - 1))) { /* not a power of 2 */
                /* what do we do here?... we need to change either
                 * f_max or deltaF so that chirplen is a power of 2
                 * so that the FFT can be done, so choose to change f_max */
                /* round chirplen up to next power of 2 */
                frexp(chirplen, &chirplen_exp);
                chirplen = ldexp(1.0, chirplen_exp);
                XLAL_PRINT_WARNING("f_max/deltaF = %g/%g = %g is not a power of two: changing f_max to %g", f_max, deltaF, f_max/deltaF, (chirplen / 2) * deltaF);
                f_max = (chirplen / 2) * deltaF;
            }
        }

        /* resize waveforms to the required length */
        XLALResizeREAL8TimeSeries(hplus, hplus->data->length - (size_t) chirplen, (size_t) chirplen);
        XLALResizeREAL8TimeSeries(hcross, hcross->data->length - (size_t) chirplen, (size_t) chirplen);

        /* put the waveform in the frequency domain */
        /* (the units will correct themselves) */
        *hptilde = XLALCreateCOMPLEX16FrequencySeries("FD H_PLUS", &hplus->epoch, 0.0, deltaF, &lalDimensionlessUnit, (size_t) chirplen / 2 + 1);
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD H_CROSS", &hcross->epoch, 0.0, deltaF, &lalDimensionlessUnit, (size_t) chirplen / 2 + 1);
        plan = XLALCreateForwardREAL8FFTPlan((size_t) chirplen, 0);
        XLALREAL8TimeFreqFFT(*hctilde, hcross, plan);
        XLALREAL8TimeFreqFFT(*hptilde, hplus, plan);

        /* clean up */
        XLALDestroyREAL8FFTPlan(plan);
        XLALDestroyREAL8TimeSeries(hcross);
        XLALDestroyREAL8TimeSeries(hplus);

    } else /* error: neither a FD nor a TD approximant */
        XLAL_ERROR(XLAL_EINVAL, "Invalid approximant");

    return 0;
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
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralChooseTDWaveform");

    return XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT, m1,
        m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i, lambda1, lambda2,
        waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
}

/** @} */

/**
 * @name General Waveform Switching Mode Generation Routines
 * @{
 */

/**
 * Interface to compute a set of -2 spin-weighted spherical harmonic modes
 * for a binary inspiral of any available amplitude and phase PN order.
 * The phasing is computed with any of the TaylorT1, T2, T3, T4 methods.
 *
 * It can also return the (2,2), (2,1), (3,3), (4,4), (5,5) modes of the EOBNRv2
 * model. Note that EOBNRv2 will ignore ampO, phaseO, lmax and f_ref arguments.
 */
SphHarmTimeSeries *XLALSimInspiralChooseTDModes(
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
    int lmax,                                   /**< generate all modes with l <= lmax */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 v0 = 1.;
    SphHarmTimeSeries *hlm = NULL;

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
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT1PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT2PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT3PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT4PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case EOBNRv2:
        case EOBNRv2HM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimIMREOBNRv2Modes(phiRef, deltaT, m1, m2, f_min, r);
            // EOB driver only outputs modes with m>0, add m<0 modes by symmetry
            size_t l, j;
            int m;
            for( l=2; l <= XLALSphHarmTimeSeriesGetMaxL( hlm ); l++ ) {
                for( m=-l; m<0; m++){
                    COMPLEX16TimeSeries* inmode = XLALSphHarmTimeSeriesGetMode(
                            hlm, l, -m );
                    if( !inmode ) continue;
                    COMPLEX16TimeSeries* tmpmode = XLALCutCOMPLEX16TimeSeries(
                            inmode, 0, inmode->data->length );
                    for( j=0; j < tmpmode->data->length; j++ ){
                        tmpmode->data->data[j] = cpow(-1, l)
                            * conj( tmpmode->data->data[j] );
                    }
                    hlm = XLALSphHarmTimeSeriesAddMode( hlm, tmpmode, l, m );
                    XLALDestroyCOMPLEX16TimeSeries( tmpmode );
                }
            }
            break;

        default:
            XLALPrintError("Cannot generate modes for this approximant\n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlm;
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
 * @param phiRef      Reference orbital phase (rad)
 * @param deltaT      Sampling interval (s)
 * @param m1          Mass of companion 1 (kg)
 * @param m2          Mass of companion 2 (kg)
 * @param f_min       Starting GW frequency (Hz)
 * @param f_ref       Reference GW frequency (Hz)
 * @param r           Distance of source (m)
 * @param lambda1     (tidal deformability of mass 1) / m1^5 (dimensionless)
 * @param lambda2     (tidal deformability of mass 2) / m2^5 (dimensionless)
 * @param waveFlags   Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags
 * @param nonGRparams Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms
 * @param amplitudeO  Twice post-Newtonian amplitude order
 * @param phaseO      Twice post-Newtonian order
 * @param lmax        Generate all modes with l <= lmax
 * @param approximant Post-Newtonian approximant to use for waveform production
 * @return Linked list of SphHarmTimeSeries modes.
 */
SphHarmTimeSeries *XLALSimInspiralModesTD(REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, int lmax, Approximant approximant)
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

    modes = XLALSimInspiralChooseTDModes(phiRef, deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, lmax, approximant);
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
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
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
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT1PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT2PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT3PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT4PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case EOBNRv2:
        case EOBNRv2HM:
            ts = XLALSimIMREOBNRv2Modes(phiRef, deltaT, m1, m2, f_min, r);
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
            /* case LAL_PNORDER_ONE: */
            case 2:
                hp1 = cos(2.*phi)*(19./6. + 3./2.*ci2 - ci4/3.
                        + eta*(-19./6. + ci2*11./6. + ci4))
                        - cos(4.*phi) * (4./3.*si2*(1. + ci2)*(1. - 3.*eta));
                hc1 = sin(2.*phi)*ci*(17./3. - ci2*4./3.
                        + eta*(-13./3. + 4.*ci2))
                        + sin(4.*phi)*ci*si2*(-8./3.*(1. - 3.*eta));
            /*case LAL_PNORDER_HALF:*/
            case 1:
                hp05 = - si*dm*(cos(phi)*(5./8. + ci2/8.)
                        - cos(3.*phi)*(9./8. + 9.*ci2/8.));
                hc05 = si*ci*dm*(-sin(phi)*3./4. + sin(3.*phi)*9./4.);
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
    REAL8 psi /**< polarization angle (rad) */
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
        // This function adds hlm(t) * Y_lm(incl,psi) to (h+ - i hx)(t)
        ret = XLALSimAddMode(*hp, *hc, ts->mode, iota, psi, ts->l, ts->m, 0);
        if( ret != XLAL_SUCCESS ) XLAL_ERROR(XLAL_EFUNC);
        ts = ts->next;
    }

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
	REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
	INT4 ampO	 	  /**< twice amp. post-Newtonian order */
	)
{
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z, lnhx, lnhy, lnhz;
    REAL8 e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, lx, ly, lz;
    REAL8 nx2, ny2, nz2, nz3, lx2, ly2, lz2, lz3;
    REAL8 hplus0, hcross0, hplus05, hcross05, hplus1, hcross1;
    REAL8 hplus15, hcross15, hplusSpin1, hcrossSpin1;
    REAL8 hplusSpin15, hcrossSpin15, hplusTail15, hcrossTail15;
    REAL8 M, eta, dm, phi, v, v2, dist, ampfac, logfac = 0.;
    INT4 idx, len;

    /* Macros to check time series vectors */
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
                logfac = log(v/v0);
                hplusTail15 = 2*((lx2 - ly2 - nx2 + ny2)*LAL_PI
                        + 12*(lx*nx - ly*ny)*logfac);
                hcrossTail15 = 4*((lx*ly - nx*ny)*LAL_PI
                        + 6*(ly*nx + lx*ny)*logfac);

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

            /*case LAL_PNORDER_HALF:*/
            case 1:
                /* 0.5PN non-spinning amp. corrections */
                hplus05 = (dm*(-2*lx2*lz + 2*ly2*lz + lz*(nx2 - ny2)
                        + 6*lx*nx*nz - 6*ly*ny*nz))/2.;
                hcross05 = dm*(-2*lx*ly*lz + lz*nx*ny
					+ 3*ly*nx*nz + 3*lx*ny*nz);

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
        case 2:
          fact1 = v4*0.25;

          *hplus += fact1*(((a2y - a1y)*e1x - (a1x - a2x)*e1y) + I*((a2y
          - a1y)*e2x - (a1x - a2x)*e2y));
          *hcross += fact1*(((a1x - a2x)*e1x - (a1y - a2y)*e1y) + I*((a1x
          - a2x)*e2x - (a1y - a2y)*e2y));
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
        case SpinTaylorT2:
        case SpinTaylorT4:
        case SpinTaylorT1:
        case IMRPhenomB:
        case PhenSpinTaylor:
        case IMRPhenomC:
	case IMRPhenomD:
	case IMRPhenomPv2:
        case PhenSpinTaylorRD:
        case SEOBNRv1:
        case SpinDominatedWf:
        case SEOBNRv2:
        case SEOBNRv2_opt:
        case SEOBNRv3:
        case SEOBNRv4:
        case SEOBNRv4_opt:
        case NR_hdf5:
            return 1;

        default:
            return 0;
    }
}

/**
 * Checks whether the given approximant is implemented in lalsimulation's XLALSimInspiralChooseFDWaveform().
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
        case IMRPhenomP:
        case IMRPhenomPv2:
        case EOBNRv2_ROM:
        case EOBNRv2HM_ROM:
        case SEOBNRv1_ROM_EffectiveSpin:
        case SEOBNRv1_ROM_DoubleSpin:
        case SEOBNRv2_ROM_EffectiveSpin:
        case SEOBNRv2_ROM_DoubleSpin:
        case SEOBNRv2_ROM_DoubleSpin_HI:
        //case TaylorR2F4:
        case TaylorF2:
	case EccentricFD:
        case SpinTaylorF2:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
        case SpinTaylorT4Fourier:
        case SpinTaylorT2Fourier:
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
    found_approximant = DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, lalSimulationApproximantNames);
    found_order       = DELETE_SUBSTRING_IN_LIST_FROM_STRING(string, lalSimulationPNOrderNames);
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
    if (XLALSimInspiralDecomposeWaveformString(&approximant, NULL, NULL, waveform) < 0)
        XLAL_ERROR(XLAL_EFUNC);
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
    case SpinTaylorT4:
    case SpinTaylorT2:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case SpinTaylorT3:
    case IMRPhenomP:
    case IMRPhenomPv2:
    case SpinTaylorT2Fourier:
    case SpinTaylorT4Fourier:
    case SpinDominatedWf:
    case SEOBNRv3:
    case NR_hdf5:
      spin_support=LAL_SIM_INSPIRAL_PRECESSINGSPIN;
      break;
    case SpinTaylorF2:
    case FindChirpPTF:
      spin_support=LAL_SIM_INSPIRAL_SINGLESPIN;
      break;
    case TaylorF2:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomB:
    case IMRPhenomC:
    case IMRPhenomD:
    case SEOBNRv1:
    case SEOBNRv2:
    case SEOBNRv4:
    case SEOBNRv2_opt:
    case SEOBNRv4_opt:
    case SEOBNRv1_ROM_EffectiveSpin:
    case SEOBNRv1_ROM_DoubleSpin:
    case SEOBNRv2_ROM_EffectiveSpin:
    case SEOBNRv2_ROM_DoubleSpin:
    case SEOBNRv2_ROM_DoubleSpin_HI:
    case TaylorR2F4:
    case IMRPhenomFB:
    case FindChirpSP:
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
      spin_support=LAL_SIM_INSPIRAL_SPINLESS;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimuation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }

    return spin_support;

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
    case SpinTaylorT2:
    case SpinTaylorT3:
    case SpinTaylorT4:
    case SpinTaylorT5:
    case SpinTaylorFrameless:
    case SpinTaylor:
    case SpinQuadTaylor:
    case FindChirpSP:
    case FindChirpPTF:
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
    case SEOBNRv1:
    case SEOBNRv2:
    case SEOBNRv2_opt:
    case SEOBNRv4:
    case SEOBNRv4_opt:
    case SEOBNRv1_ROM_EffectiveSpin:
    case SEOBNRv1_ROM_DoubleSpin:
    case SEOBNRv2_ROM_EffectiveSpin:
    case SEOBNRv2_ROM_DoubleSpin:
    case SEOBNRv2_ROM_DoubleSpin_HI:
    case IMRPhenomA:
    case IMRPhenomB:
    case IMRPhenomFA:
    case IMRPhenomFB:
    case IMRPhenomFC:
    case SpinTaylorT2Fourier:
    case SpinTaylorT4Fourier:
    case TaylorEt:
    case TaylorT4:
    case TaylorN:
    case SpinDominatedWf:
    case NumApproximants:
    case NR_hdf5:
      testGR_accept=LAL_SIM_INSPIRAL_NO_TESTGR_PARAMS;
      break;
    case TaylorF2:
    case SpinTaylorF2:
    case EccentricFD:
    case Eccentricity:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case EccentricTD:
    case IMRPhenomC:
    case IMRPhenomD:
    case IMRPhenomP:
    case IMRPhenomPv2:
      testGR_accept=LAL_SIM_INSPIRAL_TESTGR_PARAMS;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimuation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }
  return testGR_accept;
};

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

        case fSEOBNRv1Peak:
        case fSEOBNRv2Peak:
            if ( freqFunc == fSEOBNRv1Peak ) SpinAlignedEOBVersion = 1;
            if ( freqFunc == fSEOBNRv2Peak ) SpinAlignedEOBVersion = 2;
            freq = XLALSimIMRSpinAlignedEOBPeakFrequency(m1, m2, S1z, S2z,
                    SpinAlignedEOBVersion);
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
	    /* Check that spins are zero */
	    if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
	    {
		XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
		XLAL_ERROR(XLAL_EINVAL);
	    }
	case EccentricFD:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
        case TaylorF2:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
	            XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            /* Schwarzschild ISCO */
            freqFunc = fSchwarzISCO;
            break;

        /* IMR models */
        case EOBNRv2HM:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fEOBNRv2HMRD;
            break;

        case EOBNRv2:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fEOBNRv2RD;
            break;

        case SEOBNRv1:
            /* Check that the transverse spins are zero */
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
                XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fSEOBNRv1RD;
            break;

        case SEOBNRv2:
        case SEOBNRv4:
        case SEOBNRv4_opt:
        case SEOBNRv2_opt:
            /* Check that the transverse spins are zero */
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
                XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fSEOBNRv2RD;
            break;

        case IMRPhenomA:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fIMRPhenomAFinal;
            break;

        case IMRPhenomB:
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
                XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fIMRPhenomBFinal;
            break;

        case IMRPhenomC:
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            {
                XLALPrintError("Non-zero transverse spins were given, but this is a non-precessing approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
            freqFunc = fIMRPhenomCFinal;
            break;

        // FIXME: Following I don't know how to calculate */
        /* Spinning inspiral-only time domain */
        case SpinTaylorT2:
        case SpinTaylorT4:
        case SpinTaylorT1:
	case PhenSpinTaylor:
        /* Spinning with ringdown attachment */
        case PhenSpinTaylorRD:
        /* Spinning inspiral-only frequency domain */
        case SpinTaylorF2:
        /* NR waveforms */
        case NR_hdf5:
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
      if (approximant == SpinTaylorT2 || approximant == SpinTaylorT4)
          ampOrder = MAX_PRECESSING_AMP_PN_ORDER;
      else
          ampOrder = MAX_NONPRECESSING_AMP_PN_ORDER;
  }

    REAL8 fStart;
    fStart = fLow * 2./(ampOrder+2);
    return fStart;
}

/** @} */

/** @} */
