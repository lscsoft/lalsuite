/*
 *  Copyright (C) 2017 Sebastian Khan
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
 * \author Sebastian Khan
 *
 * \file
 *
 * \brief PhenomPv3HM model
 *
 * Inspiral-merger and ringdown phenomenological, frequecny domain
 * waveform model for spinning precessing binary black holes systems.
 * Models not only the dominant (l,|m|) = (2,2) modes
 * but also some of the sub-domant modes too in the co-precessing frame.
 *
 * The model has been validated against precessing NR simulations up to mass-ratio 6
 * but due to lack of precessing NR simulations above mass-ratio 6 we cannot validate it's accuracy.
 *
 * Tested Range: up to mass-ratio 6, any spin magnitude and orientation.
 * Usage Range: up to mass-ratio 20, any spin magnitude and orientation.
 */

#include <lal/Date.h>
#include <lal/Sequence.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>

#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomUtils.h"

#include "LALSimIMRPhenomPv3HM.h"

#define L_MAX_PLUS_1 5
#define PHENOM_DEFAULT_MAXF 0.5

/**
 * read in a LALDict.
 * If ModeArray in LALDict is NULL then create a ModrArray
 * with the default modes in PhenomPv3HM.
 * If ModeArray is not NULL then use the modes supplied by user.
 */
static LALDict *IMRPhenomPv3HM_setup_mode_array(
    LALDict *extraParams)
{

    /* setup ModeArray */
    if (extraParams == NULL)
        extraParams = XLALCreateDict();
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams);
    if (ModeArray == NULL)
    { /* Default behaviour */
        /* TODO: Move this into a function */
        XLAL_PRINT_INFO("Using default modes for PhenomPv3HM.\n");
        ModeArray = XLALSimInspiralCreateModeArray();
        /* Only need to define the positive m modes/
         * The negative m modes are automatically added.
         */
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 3);
        // XLALSimInspiralModeArrayPrintModes(ModeArray);
        /* Don't forget to insert ModeArray back into extraParams. */
        XLALSimInspiralWaveformParamsInsertModeArray(extraParams, ModeArray);
    }
    else
    {
        XLAL_PRINT_INFO("Using custom modes for PhenomPv3HM.\n");
    }

    XLALDestroyValue(ModeArray);
    /*TODO: Add an error check here somehow?*/

    return extraParams;
}

/**
 * Reads in a ModeArray and checks that it is valid.
 * i.e., that it contains the 2,2 mode
 * and may only contain the modes in the model
 * i.e., 22, 21, 33, 32, 44, 43
 * Only checks upto ell=8 though.
 */
static int IMRPhenomPv3HM_check_mode_array(LALValue *ModeArray)
{
    //    if 22 mode not active  -> error
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 2, 2) != 1)
    {
        XLAL_ERROR(XLAL_EFUNC, "(2,2) mode required\n");
    }

    // if no 22,21,33,32,44,43 mode and active  -> error
    // these modes are not in the model
    for (INT4 ell = 2; ell <= 8; ell++)
    {
        for (INT4 mm = -ell; mm < ell + 1; mm++)
            {
                if (ell == 2 && mm == 2){
                    continue;
                }
                else if (ell == 2 && mm == 1)
                {
                    continue;
                }
                else if (ell == 3 && mm == 3)
                {
                    continue;
                }
                else if (ell == 3 && mm == 2)
                {
                    continue;
                }
                else if (ell == 4 && mm == 4)
                {
                    continue;
                }
                else if (ell == 4 && mm == 3)
                {
                    continue;
                }

                if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, mm) == 1)
                {
                    XLAL_ERROR(XLAL_EFUNC, "(%i,%i) mode in ModeArray but model does not include this!\n", ell, mm);
                }
            }
    }
    return XLAL_SUCCESS;
}

/**
 * Precomputes useful quantities and populates the
 * PhenomPv3HMStorage and sysq (for precession angles) structs.
 */
UNUSED static int init_PhenomPv3HM_Storage(
    PhenomPv3HMStorage *p,   /**< [out] PhenomPv3Storage struct */
    sysq *pAngles,           /**< [out] precession angle pre-computations struct */
    REAL8 m1_SI,             /**< mass of primary in SI (kg) */
    REAL8 m2_SI,             /**< mass of secondary in SI (kg) */
    REAL8 S1x,               /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S1y,               /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S1z,               /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S2x,               /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 S2y,               /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 S2z,               /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 distance,    /**< distance of source (m) */
    const REAL8 inclination, /**< inclination of source (rad) */
    const REAL8 phiRef,      /**< reference orbital phase (rad) */
    const REAL8 deltaF,      /**< Sampling frequency (Hz) */
    const REAL8 f_min,       /**< Starting GW frequency (Hz) */
    const REAL8 f_max,       /**< End frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 f_ref        /**< Reference GW frequency (Hz) */
)
{
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
    XLAL_CHECK(0 != pAngles, XLAL_EFAULT, "pAngles is NULL");

    // We check if the systems is precessing because we skip angle
    // computation if this is the case.
    p->PRECESSING = 0;
    if (S1x == 0. && S1y == 0. && S2x == 0. && S2y == 0.)
    {
        p->PRECESSING = 1; // This means the system is not precessing
    }

    /* input parameters */
    p->m1_SI = m1_SI;
    p->m2_SI = m2_SI;
    p->chi1x = S1x;
    p->chi1y = S1y;
    p->chi1z = S1z;
    p->chi2x = S2x;
    p->chi2y = S2y;
    p->chi2z = S2z;
    p->distance_SI = distance;
    p->phiRef = phiRef;
    p->deltaF = deltaF;
    p->f_min = f_min;
    p->f_max = f_max;
    p->f_ref = f_ref;

    int retcode = 0;
    retcode = PhenomInternal_PrecessingSpinEnforcePrimaryIsm1(
        &(p->m1_SI),
        &(p->m2_SI),
        &(p->chi1x),
        &(p->chi1y),
        &(p->chi1z),
        &(p->chi2x),
        &(p->chi2y),
        &(p->chi2z));
    XLAL_CHECK(
        XLAL_SUCCESS == retcode,
        XLAL_EFUNC,
        "PhenomInternal_AlignedSpinEnforcePrimaryIsm1 failed");

    /* derived parameters */
    p->m1_Msun = m1_SI / LAL_MSUN_SI;
    p->m2_Msun = m2_SI / LAL_MSUN_SI;
    p->Mtot_SI = p->m1_SI + p->m2_SI;
    p->Mtot_Msun = p->m1_Msun + p->m2_Msun;

    p->eta = p->m1_Msun * p->m2_Msun / (p->Mtot_Msun * p->Mtot_Msun);
    p->q = p->m1_Msun / p->m2_Msun; /* with m1>=m2 so q>=1 */

    /* check for rounding errors */
    if (p->eta > 0.25 || p->q < 1.0)
    {
        PhenomInternal_nudge(&(p->eta), 0.25, 1e-6);
        PhenomInternal_nudge(&(p->q), 1.0, 1e-6);
    }

    p->Msec = p->Mtot_Msun * LAL_MTSUN_SI; /* Total mass in seconds */

    p->amp0 = (p->Mtot_Msun) * LAL_MRSUN_SI * (p->Mtot_Msun) * LAL_MTSUN_SI / (p->distance_SI);

    /* chi1_l == chi1z, chi2_l == chi2z so we don't need to save them as they are duplicates */
    REAL8 chi1_l, chi2_l;

    /* rotate from LAL to PhenomP frame */
    int errcode = XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(
        &chi1_l, &chi2_l, &(p->chip), &(p->thetaJN), &(p->alpha0), &(p->phi_aligned), &(p->zeta_polariz),
        p->m1_SI, p->m2_SI, p->f_ref, p->phiRef, inclination,
        p->chi1x, p->chi1y, p->chi1z,
        p->chi2x, p->chi2y, p->chi2z, IMRPhenomPv3_V);
    XLAL_CHECK(errcode == XLAL_SUCCESS, XLAL_EFUNC, "XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame failed");

    p->inclination = p->thetaJN;

    // printf("(p->zeta_polariz) = %.16f\n", (p->zeta_polariz));

    /* compute spins in polar coordinates */
    PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(&(p->chi1_theta), &(p->chi1_phi), &(p->chi1_mag), p->chi1x, p->chi1y, p->chi1z);
    PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(&(p->chi2_theta), &(p->chi2_phi), &(p->chi2_mag), p->chi2x, p->chi2y, p->chi2z);

    if (p->PRECESSING != 1) // precessing case. compute angles
    {
        /* Initialize precession angles */
        /* evaluating the angles at the reference frequency */
        p->f_ref_Orb_Hz = 0.5 * p->f_ref; /* factor of 0.5 to go from GW to Orbital frequency */
        /* precompute everything needed to compute precession angles from LALSimInspiralFDPrecAngles.c */
        /* note that the reference frequency that you pass into InitializeSystem is the GW frequency */

        /* ExpansionOrder specifies how many terms in the PN expansion of the precession angles to use.
        * In PhenomP3 we set this to 5, i.e. all but the highest order terms.
        * */
        int ExpansionOrder = 5;
        errcode = InitializeSystem(pAngles,
                                p->m1_SI, p->m2_SI,
                                LHAT_COS_THETA, LHAT_PHI,
                                cos(p->chi1_theta), p->chi1_phi, p->chi1_mag,
                                cos(p->chi2_theta), p->chi2_phi, p->chi2_mag,
                                p->f_ref, ExpansionOrder);
        XLAL_CHECK(errcode == XLAL_SUCCESS, XLAL_EFUNC, "InitializeSystem failed");
    }

    return XLAL_SUCCESS;
}

/**
 * This version doesn't construct precessing hlm modes but instead constructs hplus, hcross directly.
 */
int XLALSimIMRPhenomPv3HMGetHplusHcross(
    UNUSED COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
    UNUSED COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
    UNUSED REAL8Sequence *freqs,               /**< frequency sequency in Hz */
    UNUSED REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    UNUSED REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    UNUSED REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    UNUSED const REAL8 distance,               /**< distance of source (m) */
    UNUSED const REAL8 inclination,            /**< inclination of source (rad) */
    UNUSED const REAL8 phiRef,                 /**< reference orbital phase (rad) */
    UNUSED const REAL8 deltaF,                 /**< Sampling frequency (Hz). To use arbitrary frequency points set deltaF <= 0. */
    UNUSED REAL8 f_ref,                        /**< Reference frequency */
    UNUSED LALDict *extraParams                /**< LALDict struct */
)
{

    XLAL_CHECK(f_ref > 0, XLAL_EDOM, "f_ref must be greater than zero.\n");

    /* Use an auxiliar laldict to not overwrite the input argument */
    LALDict *extraParams_aux;

    /* setup mode array */
    if (extraParams == NULL)
    {
        extraParams_aux = XLALCreateDict();
    }
    else{
        extraParams_aux = XLALDictDuplicate(extraParams);
    }
    extraParams_aux = IMRPhenomPv3HM_setup_mode_array(extraParams_aux);
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams_aux);
    int rcode = IMRPhenomPv3HM_check_mode_array(ModeArray);
    XLAL_CHECK(XLAL_SUCCESS == rcode, rcode, "IMRPhenomPv3HM_check_mode_array failed");

    UNUSED const REAL8 Mtot_Msun = (m1_SI + m2_SI) / LAL_MSUN_SI;

    REAL8 f_min=0.;
    REAL8 f_max=0.;
    size_t npts=0;
    REAL8Sequence *freqs_seq=NULL;

    LIGOTimeGPS tC = LIGOTIMEGPSZERO; // = {0, 0}

    if ((freqs->length == 2) && (deltaF > 0.))
    {
        // uniform frequencies
        f_min = freqs->data[0];
        f_max = freqs->data[1];

        if (f_max == 0.)
        {
            f_max = XLALSimPhenomUtilsMftoHz(PHENOM_DEFAULT_MAXF, Mtot_Msun);
        }
        npts = PhenomInternal_NextPow2(f_max / deltaF) + 1;

        freqs_seq = XLALCreateREAL8Sequence(npts);
        for (size_t j = 0; j < npts; j++)
        {
            freqs_seq->data[j] = j * deltaF;
        }
        /* shift zero frequency to avoid evaluating angles at f=0.*/
        freqs_seq->data[0] = 1e-13;
        /* coalesce at t=0 */
        /* Shift by overall length in time */
        XLAL_CHECK(
            XLALGPSAdd(&tC, -1. / deltaF),
            XLAL_EFUNC,
            "Failed to shift coalescence time to t=0,\
tried to apply shift of -1.0/deltaF with deltaF=%g.",
            deltaF);

    }
    else if ((freqs->length != 2) && (deltaF <= 0.))
    { /* else if 2. i.e. not uniformly spaced then we don't shift. */
        // custom frequencies
        f_min = freqs->data[0];
        f_max = freqs->data[freqs->length - 1]; /* Last element */
        npts = freqs->length;
        freqs_seq = XLALCreateREAL8Sequence(npts);
        for (size_t j = 0; j < npts; j++)
        {
            freqs_seq->data[j] = freqs->data[j];
        }

    } else {
        XLAL_ERROR(XLAL_EFUNC, "cannot interpret frequency bounds!\n");
    }

    /* Store useful variables and compute derived and frequency independent variables */
    PhenomPv3HMStorage *pv3HM;
    pv3HM = (PhenomPv3HMStorage *)XLALMalloc(sizeof(PhenomPv3HMStorage));

    /* Struct that stores the precession angle variables */
    sysq *pAngles;
    pAngles = (sysq *)XLALMalloc(sizeof(sysq));

    int retcode = init_PhenomPv3HM_Storage(pv3HM, pAngles,
                                           m1_SI, m2_SI,
                                           chi1x, chi1y, chi1z,
                                           chi2x, chi2y, chi2z,
                                           distance, inclination, phiRef,
                                           deltaF, f_min, f_max, f_ref);
    XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "init_PhenomPv3HM_Storage failed");



    SphHarmFrequencySeries **hlmsD = XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlmsD = NULL;
    /* XLALSimIMRPhenomHMGethlmModes currently takes only the positive m modes */
    retcode = XLALSimIMRPhenomHMGethlmModes(
        hlmsD,
        freqs,
        m1_SI,
        m2_SI,
        chi1x,
        chi1y,
        chi1z,
        chi2x,
        chi2y,
        chi2z,
        pv3HM->phiRef,
        deltaF,
        f_ref,
        extraParams_aux);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "XLALSimIMRPhenomHMGethlmModes failed");

    /* Allocate hptilde and hctilde */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);
    if (!(hptilde))
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
    XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);

    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);
    if (!(hctilde))
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
    XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

    /* Frequency domain amplitude pre-factor */
    const REAL8 amp0 = XLALSimPhenomUtilsFDamp0(Mtot_Msun, distance);

    /* logic for looping over modes
    if 22 mode active:
        # compute Ylms
        for i in frequencies:
            # compute wigner-D

    if 21 mode active:
        # compute Ylms
        for i in frequencies:
            # compute wigner-D

    if 33 mode active:
        # compute Ylms
        for i in frequencies:
            # compute wigner-D

    etc...
    */

    UINT4 ell=0;
    INT4 mprime=0;

    int ret = 0;
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 2, 2) == 1)
    {
        ell=2;
        mprime=2;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);

    }
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 2, 1) == 1)
    {
        ell=2;
        mprime=1;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);
    }
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 3) == 1)
    {
        ell = 3;
        mprime = 3;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);
    }
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 2) == 1)
    {
        ell = 3;
        mprime = 2;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);
    }
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, 4) == 1)
    {
        ell = 4;
        mprime = 4;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);
    }
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, 3) == 1)
    {
        ell = 4;
        mprime = 3;

        ret = IMRPhenomPv3HM_Compute_Mode(hptilde, hctilde, ell, mprime, Mtot_Msun, pv3HM, hlmsD, pAngles, freqs_seq);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "IMRPhenomPv3HM_Compute_Mode failed for the %i, %i mode\n", ell, mprime);
    }

    // LALFree(hlmD);
    XLALDestroyREAL8Sequence(freqs_seq);
    XLALDestroyValue(ModeArray);
    XLALDestroySphHarmFrequencySeries(*hlmsD);
    XLALFree(hlmsD);

    COMPLEX16 PhPpolp, PhPpolc; /* for polarisation application */

    REAL8 cos_2zeta, sin_2zeta;
    /* apply polarisation angle - Moved this from LALSimInspiral to indide here for PhenomPv3 */
    cos_2zeta = cos(2.0 * pv3HM->zeta_polariz);
    sin_2zeta = sin(2.0 * pv3HM->zeta_polariz);

    for (UINT4 i = 0; i < (*hptilde)->data->length; i++)
    { // loop over frequency points in sequence

        // also applying Frequency Domain distance scale here - amp0
        PhPpolp = amp0 * (*hptilde)->data->data[i];
        PhPpolc = amp0 * (*hctilde)->data->data[i];
        (*hptilde)->data->data[i] = cos_2zeta * PhPpolp + sin_2zeta * PhPpolc;
        (*hctilde)->data->data[i] = cos_2zeta * PhPpolc - sin_2zeta * PhPpolp;
    }

    /* free pointers */
    LALFree(pv3HM);
    LALFree(pAngles);
    XLALDestroyDict(extraParams_aux);

    return XLAL_SUCCESS;
}

/**
 * Driver routine to compute the precessing inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomPv3 in the frequency domain.
 *
 * Reference:
 * - Hannam et al., arXiv:1308.3271 [gr-qc]
 * - Bohe et al., PhenomPv2 technical document LIGO-T1500602
 * - Chatziioannou et al., arXiv 1703.03967 [gr-qc]
 *
 * This function can be used for equally-spaced frequency series.
 * For unequal spacing, use XLALSimIMRPhenomPv3FrequencySequence instead.
 *
 * This function calls XLALSimIMRPhenomPv3HMGetHplusHcross with just
 * the l=m=2 mode
 *
 */
int XLALSimIMRPhenomPv3(
    COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
    COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
    REAL8Sequence *freqs,               /**< frequency sequency in Hz */
    REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    REAL8 S1x,                          /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S1y,                          /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S1z,                          /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 S2x,                          /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 S2y,                          /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 S2z,                          /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 distance,               /**< distance of source (m) */
    const REAL8 inclination,            /**< inclination of source (rad) */
    const REAL8 phiRef,                 /**< reference orbital phase (rad) */
    const REAL8 deltaF,                 /**< Sampling frequency (Hz) */
    const REAL8 f_ref,                  /**< Reference frequency */
    LALDict *extraParams)               /**<linked list containing the extra testing GR parameters */
{

    /* Use an auxiliar laldict to not overwrite the input argument */
    LALDict *extraParams_aux;

    /* setup mode array */
    if (extraParams == NULL)
    {
        extraParams_aux = XLALCreateDict();
    }
    else{
        extraParams_aux = XLALDictDuplicate(extraParams);
    }
    LALValue *ModeArray = XLALSimInspiralCreateModeArray();
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralWaveformParamsInsertModeArray(extraParams_aux, ModeArray);
    XLALDestroyValue(ModeArray);

    int ret = XLALSimIMRPhenomPv3HMGetHplusHcross(
        hptilde,
        hctilde,
        freqs,
        m1_SI,
        m2_SI,
        S1x,
        S1y,
        S1z,
        S2x,
        S2y,
        S2z,
        distance,
        inclination,
        phiRef,
        deltaF,
        f_ref,
        extraParams_aux);
    XLAL_CHECK(
        XLAL_SUCCESS == ret,
        XLAL_EFUNC,
        "XLALSimIMRPhenomPv3HMGetHplusHcross failed in IMRPhenomPv3");

    XLALDestroyDict(extraParams_aux);

    return XLAL_SUCCESS;
}

/**
 * This is an internal function that returns the precession angles
 * at a single frequency
 */
static int IMRPhenomPv3HM_Compute_a_b_e(REAL8 *alpha, REAL8 *beta, REAL8 *mprime_epsilon, REAL8 fHz, INT4 mprime, const REAL8 twopi_Msec, PhenomPv3HMStorage *pv3HM, sysq *pAngles)
{

    vector angles;
    REAL8 f_mprime = 0.;
    REAL8 xi = 0.;

    f_mprime = fHz / mprime;
    xi = pow(f_mprime * twopi_Msec, pAngles->onethird);
    angles = compute_phiz_zeta_costhetaL3PN(xi, pAngles);

    *alpha = angles.x + pv3HM->alpha0;
    REAL8 epsilon = angles.y;

    *mprime_epsilon = mprime * epsilon;

    /* angles.z can sometimes nan just above 1.
    * The following nudge seems to be a good fix.
    */
    PhenomInternal_nudge(&(angles.z), 1.0, 1e-6); //This could even go from 1e-6 to 1e-15 - then would have to also change in Pv3 code. Should just fix the problem in the angles code.
    *beta = acos(angles.z);

    return XLAL_SUCCESS;
}

/**
 * This is an internal function computes terms
 * required to compute hptilde and hctilde
 */
static int IMRPhenomPv3HM_wigner_loop(COMPLEX16 *Term1, COMPLEX16 *Term2, INT4 ell, INT4 mprime, IMRPhenomPv3HMYlmStruct *ylms, IMRPhenomPv3HMAlphaStruct *als, IMRPhenomPv3HMWignderStruct *wigs)
{
    INT4 minus1l=0; /* (-1)^ell */
    if (ell % 2)
        minus1l = -1;
    else
        minus1l = 1;

    COMPLEX16 Term1_sum = 0.;
    COMPLEX16 Term2_sum = 0.;

    if (ell == 2 && mprime == 2)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha2[mm + ell] * wigs->d22[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha2[-mm + ell] * wigs->d22[1][mm + ell];

            Term1_sum += ylms->Ylm2[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm2[1][mm + ell] * WigDmConj;
        }
    }
    else if (ell == 2 && mprime == 1)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha2[mm + ell] * wigs->d21[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha2[-mm + ell] * wigs->d21[1][mm + ell];

            Term1_sum += ylms->Ylm2[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm2[1][mm + ell] * WigDmConj;
        }
    }
    else if (ell == 3 && mprime == 3)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha3[mm + ell] * wigs->d33[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha3[-mm + ell] * wigs->d33[1][mm + ell];

            Term1_sum += ylms->Ylm3[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm3[1][mm + ell] * WigDmConj;
        }
    }
    else if (ell == 3 && mprime == 2)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha3[mm + ell] * wigs->d32[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha3[-mm + ell] * wigs->d32[1][mm + ell];

            Term1_sum += ylms->Ylm3[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm3[1][mm + ell] * WigDmConj;
        }
    }
    else if (ell == 4 && mprime == 4)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha4[mm + ell] * wigs->d44[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha4[-mm + ell] * wigs->d44[1][mm + ell];

            Term1_sum += ylms->Ylm4[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm4[1][mm + ell] * WigDmConj;
        }
    }
    else if (ell == 4 && mprime == 3)
    {
        for (int mm = -ell; mm <= ell; mm++)
        {
            COMPLEX16 WigD = als->alpha4[mm + ell] * wigs->d43[0][mm + ell];
            COMPLEX16 WigDmConj = als->alpha4[-mm + ell] * wigs->d43[1][mm + ell];

            Term1_sum += ylms->Ylm4[0][mm + ell] * WigD;
            Term2_sum += minus1l * ylms->Ylm4[1][mm + ell] * WigDmConj;
        }
    }
    else
    {
        XLAL_ERROR(XLAL_EFUNC, "mode %i, %i not available.\n", ell, mprime);
    }

    *Term1 = Term1_sum;
    *Term2 = Term2_sum;

    return XLAL_SUCCESS;
}

/**
 * This is an internal function that returns hptilde and hctilde
 * for a single mode in the inertial frame
 */
static int IMRPhenomPv3HM_Compute_Mode(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    UINT4 ell,
    INT4 mprime,
    const REAL8 Mtot_Msun,
    PhenomPv3HMStorage *pv3HM,
    SphHarmFrequencySeries **hlmsD,
    sysq *pAngles,
    REAL8Sequence *freqs_seq)
{

    if (pv3HM->PRECESSING == 1) // non-precessing case. skip angles
    {
        INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */
        COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlmsD, ell, mprime);
        if (!(hlm))
            XLAL_ERROR(XLAL_EFUNC);

        /* We test for hypothetical m=0 modes */
        if (mprime == 0)
        {
            sym = 0;
        }
        else
        {
            sym = 1;
        }
        PhenomInternal_IMRPhenomHMFDAddMode(*hptilde, *hctilde, hlm, pv3HM->inclination, 0., ell, mprime, sym);
    } else { // precessing case. compute angles and do the twist

        const REAL8 Msec = Mtot_Msun * LAL_MTSUN_SI; /* Total mass in seconds */
        const REAL8 twopi_Msec = LAL_TWOPI * Msec;
        REAL8 fHz = 0.;

        COMPLEX16 half_amp_eps;

        COMPLEX16 Term1_sum = 0;
        COMPLEX16 Term2_sum = 0;

        UNUSED INT4 minus1l = 0; /* (-1)^ell */
        int ret_abe;
        int retloop;

        REAL8 alpha = 0.;
        REAL8 beta = 0.;
        REAL8 mprime_epsilon = 0.;

        // compute Ylms
        IMRPhenomPv3HMYlmStruct *ylms = XLALSimIMRPhenomPv3HMComputeYlmElements(pv3HM->inclination, 0, ell);
        // get non-prec mode
        COMPLEX16FrequencySeries *hlmD = XLALSphHarmFrequencySeriesGetMode(*hlmsD, ell, mprime);
        if (!(hlmD))
            XLAL_ERROR(XLAL_EFUNC);

        IMRPhenomPv3HMAlphaStruct *als = XLALMalloc(sizeof(IMRPhenomPv3HMAlphaStruct));
        memset(als, 0, sizeof(IMRPhenomPv3HMAlphaStruct));

        IMRPhenomPv3HMWignderStruct *wigs = XLALMalloc(sizeof(IMRPhenomPv3HMWignderStruct));
        memset(wigs, 0, sizeof(IMRPhenomPv3HMWignderStruct));

        int ret_als;
        int ret_wigs;

        // frequency loop
        for (size_t j = 0; j < freqs_seq->length; j++)
        {
            fHz = freqs_seq->data[j]; //for the angles
            // compute alpha, beta and mprime*epsilon
            ret_abe = IMRPhenomPv3HM_Compute_a_b_e(&alpha, &beta, &mprime_epsilon, fHz, mprime, twopi_Msec, pv3HM, pAngles);
            XLAL_CHECK(
                XLAL_SUCCESS == ret_abe,
                XLAL_EFUNC,
                "IMRPhenomPv3HM_Compute_a_b_e failed");
            /* Precompute wigner-d elements */
            ret_wigs = XLALSimIMRPhenomPv3HMComputeWignerdElements(&wigs, ell, mprime, -beta);
            XLAL_CHECK(
                XLAL_SUCCESS == ret_wigs,
                XLAL_EFUNC,
                "XLALSimIMRPhenomPv3HMComputeWignerdElements failed");
            /* Precompute powers of e^{i m alpha} */
            ret_als = XLALSimIMRPhenomPv3HMComputeAlphaElements(&als, ell, alpha);
            XLAL_CHECK(
                XLAL_SUCCESS == ret_als,
                XLAL_EFUNC,
                "XLALSimIMRPhenomPv3HMComputeAlphaElements failed");
            retloop = IMRPhenomPv3HM_wigner_loop(&Term1_sum, &Term2_sum, ell, mprime, ylms, als, wigs);
            XLAL_CHECK(
                XLAL_SUCCESS == retloop,
                XLAL_EFUNC,
                "IMRPhenomPv3HM_wigner_loop failed");

            COMPLEX16 hlmD_j = (hlmD->data->data[j]);
            half_amp_eps = 0.5 * hlmD_j * cexp(-I * mprime_epsilon);
            (*hptilde)->data->data[j] += half_amp_eps * (Term1_sum + Term2_sum);
            (*hctilde)->data->data[j] += -I * half_amp_eps * (Term1_sum - Term2_sum);
        }

        XLALFree(wigs);
        XLALFree(als);
        XLALFree(ylms); // allocated in XLALSimIMRPhenomPv3HMComputeYlmElements
    }

    return XLAL_SUCCESS;
}

/**
 * Returns frequency domain hlm's in inertial frame
 *
 */
int XLALSimIMRPhenomPv3HMModes(
    SphHarmFrequencySeries **hlms, /**< [out] SphHarmFrequencySeries struct containing inertial frame hlm modes */
    REAL8Sequence *freqs,               /**< frequency sequency in Hz */
    REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    REAL8 chi1x,                          /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y,                          /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z,                          /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x,                          /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y,                          /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z,                          /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 phiRef,                 /**< reference orbital phase (rad) */
    const REAL8 deltaF,                 /**< Sampling frequency (Hz) */
    const REAL8 f_ref,                  /**< Reference frequency */
    LALDict *extraParams)               /**<linked list containing the extra testing GR parameters */
{

    XLAL_ERROR(XLAL_EFUNC, "Function (XLALSimIMRPhenomPv3HMModes) not implemented!\n");

    XLAL_CHECK(f_ref > 0, XLAL_EDOM, "f_ref must be greater than zero.\n");

    /* Use an auxiliar laldict to not overwrite the input argument */
    LALDict *extraParams_aux;

    /* setup mode array */
    if (extraParams == NULL)
    {
        extraParams_aux = XLALCreateDict();
    }
    else{
        extraParams_aux = XLALDictDuplicate(extraParams);
    }
    extraParams_aux = IMRPhenomPv3HM_setup_mode_array(extraParams_aux);
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams_aux);
    int rcode = IMRPhenomPv3HM_check_mode_array(ModeArray);
    XLAL_CHECK(XLAL_SUCCESS == rcode, rcode, "IMRPhenomPv3HM_check_mode_array failed");

    UNUSED const REAL8 Mtot_Msun = (m1_SI + m2_SI) / LAL_MSUN_SI;

    REAL8 f_min = 0.;
    REAL8 f_max = 0.;
    size_t npts = 0;
    REAL8Sequence *freqs_seq = NULL;

    LIGOTimeGPS tC = LIGOTIMEGPSZERO; // = {0, 0}

    if ((freqs->length == 2) && (deltaF > 0.))
    {
        // uniform frequencies
        f_min = freqs->data[0];
        f_max = freqs->data[1];

        if (f_max == 0.)
        {
            f_max = XLALSimPhenomUtilsMftoHz(0.5, Mtot_Msun);
        }
        npts = PhenomInternal_NextPow2(f_max / deltaF) + 1;

        freqs_seq = XLALCreateREAL8Sequence(npts);
        for (size_t j = 0; j < npts; j++)
        {
            freqs_seq->data[j] = j * deltaF;
        }
        /* coalesce at t=0 */
        /* Shift by overall length in time */
        XLAL_CHECK(
            XLALGPSAdd(&tC, -1. / deltaF),
            XLAL_EFUNC,
            "Failed to shift coalescence time to t=0,\
tried to apply shift of -1.0/deltaF with deltaF=%g.",
            deltaF);
    }
    else if ((freqs->length != 2) && (deltaF <= 0.))
    { /* else if 2. i.e. not uniformly spaced then we don't shift. */
        // custom frequencies
        f_min = freqs->data[0];
        f_max = freqs->data[freqs->length - 1]; /* Last element */
        npts = freqs->length;
        freqs_seq = XLALCreateREAL8Sequence(npts);
        for (size_t j = 0; j < npts; j++)
        {
            freqs_seq->data[j] = freqs->data[j];
        }
    }
    else
    {
        XLAL_ERROR(XLAL_EDOM, "Input frequencies not as expected. Aborting.\n");
    }

    /* Store useful variables and compute derived and frequency independent variables */
    PhenomPv3HMStorage *pv3HM;
    pv3HM = (PhenomPv3HMStorage *)XLALMalloc(sizeof(PhenomPv3HMStorage));

    /* Struct that stores the precession angle variables */
    sysq *pAngles;
    pAngles = (sysq *)XLALMalloc(sizeof(sysq));

    // distance and inclination are not used in the modes
    // so we just parse dummy variables
    REAL8 distance = 1.;
    REAL8 inclination = 0.;

    int retcode = init_PhenomPv3HM_Storage(pv3HM, pAngles,
                                           m1_SI, m2_SI,
                                           chi1x, chi1y, chi1z,
                                           chi2x, chi2y, chi2z,
                                           distance, inclination, phiRef,
                                           deltaF, f_min, f_max, f_ref);
    XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "init_PhenomPv3HM_Storage failed");

    // hlmsD = co-precessing frame modes
    SphHarmFrequencySeries **hlmsD = XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlmsD = NULL;
    /* XLALSimIMRPhenomHMGethlmModes currently takes only the positive m modes */
    retcode = XLALSimIMRPhenomHMGethlmModes(
        hlmsD,
        freqs,
        m1_SI,
        m2_SI,
        chi1x,
        chi1y,
        chi1z,
        chi2x,
        chi2y,
        chi2z,
        pv3HM->phiRef,
        deltaF,
        f_ref,
        extraParams_aux);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "XLALSimIMRPhenomHMGethlmModes failed");

    const REAL8 Msec = Mtot_Msun * LAL_MTSUN_SI; /* Total mass in seconds */
    const REAL8 twopi_Msec = LAL_TWOPI * Msec;
    REAL8 fHz = 0.;

    int ret_abe;
    int ret_wig_element;

    REAL8 alpha = 0.;
    REAL8 beta = 0.;
    REAL8 mprime_epsilon = 0.;
    REAL8 wig_d = 0.;

    /* loop over modes */
    /* at this point ModeArray should contain the list of modes
     * and therefore if NULL then something is wrong and abort.
     */
    if (ModeArray == NULL)
    {
        XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
    }
    for (UINT4 ell = 2; ell < L_MAX_PLUS_1; ell++)
    { // inertial frame ell modes
        for (INT4 mm = -ell; mm < (INT4)ell + 1; mm++)
        { // inertial frame mm modes
            // inertial frame mode FrequencySeries
            COMPLEX16FrequencySeries *hlm = XLALCreateCOMPLEX16FrequencySeries("hlm: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);
            memset(hlm->data->data, 0, npts * sizeof(COMPLEX16));

            for (INT4 mprime = 1; mprime < (INT4)ell + 1; mprime++)
            {
                 /* The negative m modes are automatically added. */
                 /* first check if (l,m) mode is 'activated' in the ModeArray */
                 /* if activated then generate the mode, else skip this mode. */
                 if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, mprime) != 1)
                { /* skip mode */
                    XLAL_PRINT_INFO("SKIPPING ell = %i mprime = %i\n", ell, mprime);
                    continue;
                } /* else: generate mode */
                XLAL_PRINT_INFO("generateing ell = %i mprime = %i\n", ell, mprime);

                COMPLEX16FrequencySeries *hlmD = XLALSphHarmFrequencySeriesGetMode(*hlmsD, ell, mprime);
                if (!(hlmD))
                    XLAL_ERROR(XLAL_EFUNC, "XLALSphHarmFrequencySeriesGetMode failed for (%i,%i) mode\n", ell, mprime);

                // frequency loop
                for (size_t j = 0; j < freqs_seq->length; j++)
                {

                    fHz = freqs_seq->data[j]; //for the angles
                    // compute alpha, beta and mprime*epsilon
                    ret_abe = IMRPhenomPv3HM_Compute_a_b_e(&alpha, &beta, &mprime_epsilon, fHz, mprime, twopi_Msec, pv3HM, pAngles);
                    XLAL_CHECK(
                        XLAL_SUCCESS == ret_abe,
                        XLAL_EFUNC,
                        "IMRPhenomPv3HM_Compute_a_b_e failed");

                    ret_wig_element = XLALSimPhenomUtilsPhenomPv3HMWignerdElement(&wig_d, ell, mprime, mm, -beta);
                    XLAL_CHECK(
                        XLAL_SUCCESS == ret_wig_element,
                        XLAL_EFUNC,
                        "XLALSimPhenomUtilsPhenomPv3HMWignerdElement failed");

                    COMPLEX16 hlmD_j = hlmD->data->data[j];
                    hlm->data->data[j] += hlmD_j * cexp(-I * mprime_epsilon) * cexp(I * mm * alpha) * wig_d;
                }
            }

            *hlms = XLALSphHarmFrequencySeriesAddMode(*hlms, hlm, ell, mm);

            XLALDestroyCOMPLEX16FrequencySeries(hlm);
        }
    }

    // LALFree(hlmD);
    XLALDestroyREAL8Sequence(freqs_seq);
    XLALDestroyValue(ModeArray);
    XLALDestroySphHarmFrequencySeries(*hlmsD);
    XLALFree(hlmsD);


    /* free pointers */
    LALFree(pv3HM);
    LALFree(pAngles);
    XLALDestroyDict(extraParams_aux);

    return XLAL_SUCCESS;
}
