/*
 *  Copyright (C) 2017 Sebastian Khan, Francesco Pannarale, Lionel London
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

#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>
#include <lal/Date.h>
#include <lal/Units.h>

#include "LALSimIMRPhenomHM.h"
#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomUtils.h"
#include "LALSimRingdownCW.h"
#include "LALSimIMRPhenomD_internals.c"

/*
 * Phase shift due to leading order complex amplitude
 * [L.Blancet, arXiv:1310.1528 (Sec. 9.5)]
 * "Spherical hrmonic modes for numerical relativity"
 */
/* List of phase shifts: the index is the azimuthal number m */
static const double cShift[7] = {0.0,
                                 LAL_PI_2 /* i shift */,
                                 0.0,
                                 -LAL_PI_2 /* -i shift */,
                                 LAL_PI /* 1 shift */,
                                 LAL_PI_2 /* -1 shift */,
                                 0.0};

/**
 * read in a LALDict.
 * If ModeArray in LALDict is NULL then create a ModrArray
 * with the default modes in PhenomHM.
 * If ModeArray is not NULL then use the modes supplied by user.
 */
LALDict *IMRPhenomHM_setup_mode_array(
    LALDict *extraParams)
{

    /* setup ModeArray */
    if (extraParams == NULL)
        extraParams = XLALCreateDict();
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams);
    if (ModeArray == NULL)
    { /* Default behaviour */
        XLAL_PRINT_INFO("Using default modes for PhenomHM.\n");
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
        XLAL_PRINT_INFO("Using custom modes for PhenomHM.\n");
    }

    XLALDestroyValue(ModeArray);

    return extraParams;
}

/**
 *
 */
int PhenomHM_init_useful_mf_powers(PhenomHMUsefulMfPowers *p, REAL8 number)
{
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
    XLAL_CHECK(!(number < 0), XLAL_EDOM, "number must be non-negative");

    // consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
    p->itself = number;
    p->sqrt = sqrt(number);
    p->sixth = cbrt(p->sqrt);
    p->m_sixth = 1.0 / p->sixth;
    p->third = p->sixth * p->sixth;
    p->two_thirds = p->third * p->third;
    p->four_thirds = number * p->third;
    p->five_thirds = number * p->two_thirds;
    p->two = number * number;
    p->seven_thirds = p->third * p->two;
    p->eight_thirds = p->two_thirds * p->two;
    p->m_seven_sixths = p->m_sixth / number;
    p->m_five_sixths = p->m_seven_sixths * p->third;
    p->m_sqrt = 1. / p->sqrt;

    return XLAL_SUCCESS;
}

int PhenomHM_init_useful_powers(PhenomHMUsefulPowers *p, REAL8 number)
{
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
    XLAL_CHECK(number >= 0, XLAL_EDOM, "number must be non-negative");

    // consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
    double sixth = pow(number, 1.0 / 6.0);
    p->third = sixth * sixth;
    p->two_thirds = p->third * p->third;
    p->four_thirds = number * p->third;
    p->five_thirds = p->four_thirds * p->third;
    p->two = number * number;
    p->seven_thirds = p->third * p->two;
    p->eight_thirds = p->two_thirds * p->two;
    p->inv = 1. / number;
    double m_sixth = 1.0 / sixth;
    p->m_seven_sixths = p->inv * m_sixth;
    p->m_third = m_sixth * m_sixth;
    p->m_two_thirds = p->m_third * p->m_third;
    p->m_five_thirds = p->inv * p->m_two_thirds;

    return XLAL_SUCCESS;
}

/**
 * helper function to multiple hlm with Ylm.
 * Adapted from LALSimIMREOBNRv2HMROMUtilities.c
 */
int IMRPhenomHMFDAddMode(
    COMPLEX16FrequencySeries *hptilde,
    COMPLEX16FrequencySeries *hctilde,
    COMPLEX16FrequencySeries *hlmtilde,
    REAL8 theta,
    REAL8 phi,
    INT4 l,
    INT4 m,
    INT4 sym)
{
    COMPLEX16 Y;
    UINT4 j;
    COMPLEX16 hlm; /* helper variable that contain a single point of hlmtilde */

    INT4 minus1l; /* (-1)^l */
    if (l % 2)
        minus1l = -1;
    else
        minus1l = 1;
    if (sym)
    { /* Equatorial symmetry: add in -m mode */
        Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
        COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
        COMPLEX16 factorp = 0.5 * (Y + minus1l * Ymstar);
        COMPLEX16 factorc = I * 0.5 * (Y - minus1l * Ymstar);
        for (j = 0; j < hlmtilde->data->length; ++j)
        {
            hlm = (hlmtilde->data->data[j]);
            hptilde->data->data[j] += factorp * hlm;
            hctilde->data->data[j] += factorc * hlm;
        }
    }
    else
    { /* not adding in the -m mode */
        Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
        COMPLEX16 factorp = 0.5 * Y;
        COMPLEX16 factorc = I * factorp;
        for (j = 0; j < hlmtilde->data->length; ++j)
        {
            hlm = (hlmtilde->data->data[j]);
            hptilde->data->data[j] += factorp * hlm;
            hctilde->data->data[j] += factorc * hlm;
        }
    }

    return XLAL_SUCCESS;
}

/**
 * returns the real and imag parts of the complex ringdown frequency
 * for the (l,m) mode.
 */
int IMRPhenomHMGetRingdownFrequency(
    REAL8 *fringdown,
    REAL8 *fdamp,
    UINT4 ell,
    INT4 mm,
    REAL8 finalmass,
    REAL8 finalspin)
{
    const REAL8 inv2Pi = 0.5 / LAL_PI;
    complex double ZZ;
    ZZ = SimRingdownCW_CW07102016(SimRingdownCW_KAPPA(finalspin, ell, mm), ell, mm, 0);
    const REAL8 Mf_RD_tmp = inv2Pi * creal(ZZ); /* GW ringdown frequency, converted from angular frequency */
    *fringdown = Mf_RD_tmp / finalmass;         /* scale by predicted final mass */
    /* lm mode ringdown damping time (imaginary part of ringdown), geometric units */
    const REAL8 f_DAMP_tmp = inv2Pi * cimag(ZZ); /* this is the 1./tau in the complex QNM */
    *fdamp = f_DAMP_tmp / finalmass;             /* scale by predicted final mass */

    return XLAL_SUCCESS;
}

/**
 * helper function to easily check if the
 * input frequency sequence is uniformly space
 * or a user defined set of discrete frequencies.
 */
UINT4 IMRPhenomHM_is_freq_uniform(
    REAL8Sequence *freqs,
    REAL8 deltaF)
{
    UINT4 freq_is_uniform = 0;
    if ((freqs->length == 2) && (deltaF > 0.))
    {
        freq_is_uniform = 1;
    }
    else if ((freqs->length != 2) && (deltaF <= 0.))
    {
        freq_is_uniform = 0;
    }

    return freq_is_uniform;
}

/**
 * derive frequency variables for PhenomHM based on input.
 * used to set the index on arrays where we have non-zero values.
 */
int init_IMRPhenomHMGet_FrequencyBounds_storage(
    PhenomHMFrequencyBoundsStorage *p, /**< [out] PhenomHMFrequencyBoundsStorage struct */
    REAL8Sequence *freqs,              /**< Input list of GW frequencies [Hz] */
    REAL8 Mtot,                        /**< total mass in solar masses */
    REAL8 deltaF,                      /**< frequency spacing */
    REAL8 f_ref_in                     /**< reference GW frequency */
)
{
    p->deltaF = deltaF;
    /* determine how to populate frequency sequence */
    /* if len(freqs_in) == 2 and deltaF > 0. then
     * f_min = freqs_in[0]
     * f_max = freqs_in[1]
     * else if len(freqs_in) != 2 and deltaF <= 0. then
     * user has given an arbitrary set of frequencies to evaluate the model at.
     */

    p->freq_is_uniform = IMRPhenomHM_is_freq_uniform(freqs, p->deltaF);

    if (p->freq_is_uniform == 1)
    { /* This case we use regularly spaced frequencies */
        p->f_min = freqs->data[0];
        p->f_max = freqs->data[1];

        /* If p->f_max == 0. Then we default to the ending frequency
         * for PhenomHM
         */
        if (p->f_max == 0.)
        {
            p->f_max = XLALSimPhenomUtilsMftoHz(
                PHENOMHM_DEFAULT_MF_MAX, Mtot);
        }
        /* we only need to evaluate the phase from
         * f_min to f_max with a spacing of deltaF
         */
        p->npts = PhenomInternal_NextPow2(p->f_max / p->deltaF) + 1;
        p->ind_min = (size_t)ceil(p->f_min / p->deltaF);
        p->ind_max = (size_t)ceil(p->f_max / p->deltaF);
        XLAL_CHECK((p->ind_max <= p->npts) && (p->ind_min <= p->ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=npts=%zu.", p->ind_min, p->ind_max, p->npts);
    }
    else if (p->freq_is_uniform == 0)
    { /* This case we possibly use irregularly spaced frequencies */
        /* Check that the frequencies are always increasing */
        for (UINT4 i = 0; i < freqs->length - 1; i++)
        {
            XLAL_CHECK(
                freqs->data[i] - freqs->data[i + 1] < 0.,
                XLAL_EFUNC,
                "custom frequencies must be increasing.");
        }

        XLAL_PRINT_INFO("Using custom frequency input.\n");
        p->f_min = freqs->data[0];
        p->f_max = freqs->data[freqs->length - 1]; /* Last element */

        p->npts = freqs->length;
        p->ind_min = 0;
        p->ind_max = p->npts;
    }
    else
    { /* Throw an informative error. */
        XLAL_PRINT_ERROR("Input sequence of frequencies and deltaF is not \
    compatible.\nSpecify a f_min and f_max by using a REAL8Sequence of length = 2 \
    along with a deltaF > 0.\
    \nIf you want to supply an arbitrary list of frequencies to evaluate the with \
    then supply those frequencies using a REAL8Sequence and also set deltaF <= 0.");
    }

    /* Fix default behaviour for f_ref */
    /* If f_ref = 0. then set f_ref = f_min */
    p->f_ref = f_ref_in;
    if (p->f_ref == 0.)
    {
        p->f_ref = p->f_min;
    }

    return XLAL_SUCCESS;
}

/**
 * Precompute a bunch of PhenomHM related quantities and store them filling in a
 * PhenomHMStorage variable
 */
static int init_PhenomHM_Storage(
    PhenomHMStorage *p,
    const REAL8 m1_SI,
    const REAL8 m2_SI,
    const REAL8 chi1z,
    const REAL8 chi2z,
    REAL8Sequence *freqs,
    const REAL8 deltaF,
    const REAL8 f_ref,
    const REAL8 phiRef)
{
    int retcode;
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");

    p->m1 = m1_SI / LAL_MSUN_SI;
    p->m2 = m2_SI / LAL_MSUN_SI;
    p->m1_SI = m1_SI;
    p->m2_SI = m2_SI;
    p->Mtot = p->m1 + p->m2;
    p->eta = p->m1 * p->m2 / (p->Mtot * p->Mtot);
    p->chi1z = chi1z;
    p->chi2z = chi2z;
    p->phiRef = phiRef;
    p->deltaF = deltaF;
    p->freqs = freqs;

    if (p->eta > 0.25)
        PhenomInternal_nudge(&(p->eta), 0.25, 1e-6);
    if (p->eta > 0.25 || p->eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");
    if (p->eta < MAX_ALLOWED_ETA)
        XLAL_PRINT_WARNING("Warning: The model is not calibrated for mass-ratios above 20\n");

    retcode = 0;
    retcode = PhenomInternal_AlignedSpinEnforcePrimaryIsm1(
        &(p->m1),
        &(p->m2),
        &(p->chi1z),
        &(p->chi2z));
    XLAL_CHECK(
        XLAL_SUCCESS == retcode,
        XLAL_EFUNC,
        "PhenomInternal_AlignedSpinEnforcePrimaryIsm1 failed");

    /* sanity checks on frequencies */
    PhenomHMFrequencyBoundsStorage pHMFS;
    retcode = 0;
    retcode = init_IMRPhenomHMGet_FrequencyBounds_storage(
        &pHMFS,
        p->freqs,
        p->Mtot,
        p->deltaF,
        f_ref);
    XLAL_CHECK(
        XLAL_SUCCESS == retcode,
        XLAL_EFUNC,
        "init_IMRPhenomHMGet_FrequencyBounds_storage failed");

    /* redundent storage */
    p->f_min = pHMFS.f_min;
    p->f_max = pHMFS.f_max;
    p->f_ref = pHMFS.f_ref;
    p->freq_is_uniform = pHMFS.freq_is_uniform;
    p->npts = pHMFS.npts;
    p->ind_min = pHMFS.ind_min;
    p->ind_max = pHMFS.ind_max;

    p->Mf_ref = XLALSimPhenomUtilsHztoMf(p->f_ref, p->Mtot);

    p->finmass = IMRPhenomDFinalMass(p->m1, p->m2, p->chi1z, p->chi2z);
    p->finspin = XLALSimIMRPhenomDFinalSpin(p->m1, p->m2, p->chi1z, p->chi2z); /* dimensionless final spin */
    if (p->finspin > 1.0)
        XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    /* populate the ringdown frequency array */
    /* If you want to model a new mode then you have to add it here. */
    /* (l,m) = (2,2) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[2][2],
        &p->PhenomHMfdamp[2][2],
        2, 2,
        p->finmass, p->finspin);

    /* (l,m) = (2,1) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[2][1],
        &p->PhenomHMfdamp[2][1],
        2, 1,
        p->finmass, p->finspin);

    /* (l,m) = (3,3) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[3][3],
        &p->PhenomHMfdamp[3][3],
        3, 3,
        p->finmass, p->finspin);

    /* (l,m) = (3,2) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[3][2],
        &p->PhenomHMfdamp[3][2],
        3, 2,
        p->finmass, p->finspin);

    /* (l,m) = (4,4) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[4][4],
        &p->PhenomHMfdamp[4][4],
        4, 4,
        p->finmass, p->finspin);

    /* (l,m) = (4,3) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[4][3],
        &p->PhenomHMfdamp[4][3],
        4, 3,
        p->finmass, p->finspin);

    p->Mf_RD_22 = p->PhenomHMfring[2][2];
    p->Mf_DM_22 = p->PhenomHMfdamp[2][2];

    /* (l,m) = (2,2) */
    int ell, mm;
    ell = 2;
    mm = 2;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;
    /* (l,m) = (2,1) */
    ell = 2;
    mm = 1;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;
    /* (l,m) = (3,3) */
    ell = 3;
    mm = 3;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;
    /* (l,m) = (3,2) */
    ell = 3;
    mm = 2;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;
    /* (l,m) = (4,4) */
    ell = 4;
    mm = 4;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;
    /* (l,m) = (4,3) */
    ell = 4;
    mm = 3;
    p->Rholm[ell][mm] = p->Mf_RD_22 / p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm] / p->Mf_DM_22;

    return XLAL_SUCCESS;
};

/**
 * mathematica function postfRDflm
 * domain mapping function - ringdown
 */
double IMRPhenomHMTrd(
    REAL8 Mf,
    REAL8 Mf_RD_22,
    REAL8 Mf_RD_lm,
    const INT4 AmpFlag,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM)
{
    double ans = 0.0;
    if (AmpFlag == AmpFlagTrue)
    {
        /* For amplitude */
        ans = Mf - Mf_RD_lm + Mf_RD_22; /*Used for the Amplitude as an approx fix for post merger powerlaw slope */
    }
    else
    {
        /* For phase */
        REAL8 Rholm = pHM->Rholm[ell][mm];
        ans = Rholm * Mf; /* Used for the Phase */
    }

    return ans;
}

/**
 * mathematica function Ti
 * domain mapping function - inspiral
 */
double IMRPhenomHMTi(REAL8 Mf, const INT4 mm)
{
    return 2.0 * Mf / mm;
}

/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
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
    PhenomHMStorage *pHM)
{
    REAL8 Trd = IMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, pHM);
    REAL8 Ti = IMRPhenomHMTi(fi, mm);

    //Am = ( Trd[fr]-Ti[fi] )/( fr - fi );
    *Am = (Trd - Ti) / (fr - fi);

    //Bm = Ti[fi] - fi*Am;
    *Bm = Ti - fi * (*Am);

    return XLAL_SUCCESS;
}

/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
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
    REAL8 Br)
{
    // Define function to output map params used depending on
    if (flm > fi)
    {
        if (flm > fr)
        {
            *a = Ar;
            *b = Br;
        }
        else
        {
            *a = Am;
            *b = Bm;
        }
    }
    else
    {
        *a = Ai;
        *b = Bi;
    };
    return XLAL_SUCCESS;
}

/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
int IMRPhenomHMFreqDomainMapParams(
    REAL8 *a,             /**< [Out]  */
    REAL8 *b,             /**< [Out]  */
    REAL8 *fi,            /**< [Out]  */
    REAL8 *fr,            /**< [Out]  */
    REAL8 *f1,            /**< [Out]  */
    const REAL8 flm,      /**< input waveform frequency */
    const INT4 ell,       /**< spherical harmonics ell mode */
    const INT4 mm,        /**< spherical harmonics m mode */
    PhenomHMStorage *pHM, /**< Stores quantities in order to calculate them only once */
    const int AmpFlag     /**< is ==1 then computes for amplitude, if ==0 then computes for phase */
)
{

    /*check output points are NULL*/
    XLAL_CHECK(a != NULL, XLAL_EFAULT);
    XLAL_CHECK(b != NULL, XLAL_EFAULT);
    XLAL_CHECK(fi != NULL, XLAL_EFAULT);
    XLAL_CHECK(fr != NULL, XLAL_EFAULT);
    XLAL_CHECK(f1 != NULL, XLAL_EFAULT);

    /* Account for different f1 definition between PhenomD Amplitude and Phase derivative models */
    REAL8 Mf_1_22 = 0.; /* initalise variable */
    if (AmpFlag == AmpFlagTrue)
    {
        /* For amplitude */
        Mf_1_22 = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD [amplitude model], for the 22 mode */
    }
    else
    {
        /* For phase */
        Mf_1_22 = PHI_fJoin_INS; /* inspiral joining frequency from PhenomD [phase model], for the 22 mode */
    }

    *f1 = Mf_1_22;

    REAL8 Mf_RD_22 = pHM->Mf_RD_22;
    REAL8 Mf_RD_lm = pHM->PhenomHMfring[ell][mm];

    // Define a ratio of QNM frequencies to be used for scaling various quantities
    REAL8 Rholm = pHM->Rholm[ell][mm];

    // Given experiments with the l!=m modes, it appears that the QNM scaling rather than the PN scaling may be optimal for mapping f1
    REAL8 Mf_1_lm = Mf_1_22 / Rholm;

    /* Define transition frequencies */
    *fi = Mf_1_lm;
    *fr = Mf_RD_lm;

    /*Define the slope and intercepts of the linear transformation used*/
    REAL8 Ai = 2.0 / mm;
    REAL8 Bi = 0.0;
    REAL8 Am;
    REAL8 Bm;
    IMRPhenomHMSlopeAmAndBm(&Am, &Bm, mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, pHM);

    REAL8 Ar = 1.0;
    REAL8 Br = 0.0;
    if (AmpFlag == AmpFlagTrue)
    {
        /* For amplitude */
        Br = -Mf_RD_lm + Mf_RD_22;
    }
    else
    {
        /* For phase */
        Ar = Rholm;
    }

    /* Define function to output map params used depending on */
    int ret = IMRPhenomHMMapParams(a, b, flm, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMMapParams failed in IMRPhenomHMFreqDomainMapParams (1)\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    return XLAL_SUCCESS;
}

/**
 * IMRPhenomHMFreqDomainMap
 * Input waveform frequency in Geometric units (Mflm)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double IMRPhenomHMFreqDomainMap(
    REAL8 Mflm,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM,
    const int AmpFlag)
{

    /* Mflm here has the same meaning as Mf_wf in XLALSimIMRPhenomHMFreqDomainMapHM (old deleted function). */
    REAL8 a = 0.;
    REAL8 b = 0.;
    /* Following variables not used in this funciton but are returned in IMRPhenomHMFreqDomainMapParams */
    REAL8 fi = 0.;
    REAL8 fr = 0.;
    REAL8 f1 = 0.;
    int ret = IMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, Mflm, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMap\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Mf22 = a * Mflm + b;
    return Mf22;
}

int IMRPhenomHMPhasePreComp(
    HMPhasePreComp *q,          /**< [out] HMPhasePreComp struct */
    const INT4 ell,             /**< ell spherical harmonic number */
    const INT4 mm,              /**< m spherical harmonic number */
    PhenomHMStorage *pHM,       /**< PhenomHMStorage struct */
    UNUSED LALDict *extraParams /**< LALDict strcut */
)
{
    REAL8 ai = 0.0;
    REAL8 bi = 0.0;
    REAL8 am = 0.0;
    REAL8 bm = 0.0;
    REAL8 ar = 0.0;
    REAL8 br = 0.0;
    REAL8 fi = 0.0;
    REAL8 f1 = 0.0;
    REAL8 fr = 0.0;

    const INT4 AmpFlag = 0;

    /* NOTE: As long as Mfshit isn't >= fr then the value of the shift is arbitrary. */
    const REAL8 Mfshift = 0.0001;

    int ret = IMRPhenomHMFreqDomainMapParams(&ai, &bi, &fi, &fr, &f1, Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - inspiral\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    q->ai = ai;
    q->bi = bi;

    ret = IMRPhenomHMFreqDomainMapParams(&am, &bm, &fi, &fr, &f1, fi + Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - intermediate\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    q->am = am;
    q->bm = bm;

    ret = IMRPhenomHMFreqDomainMapParams(&ar, &br, &fi, &fr, &f1, fr + Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - merger-ringdown\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ar = ar;
    q->br = br;

    q->fi = fi;
    q->fr = fr;

    double Rholm = pHM->Rholm[ell][mm];
    double Taulm = pHM->Taulm[ell][mm];

    PhenDAmpAndPhasePreComp pDPreComp;
    ret = IMRPhenomDSetupAmpAndPhaseCoefficients(
        &pDPreComp,
        pHM->m1,
        pHM->m2,
        pHM->chi1z,
        pHM->chi2z,
        Rholm,
        Taulm,
        extraParams);
    if (ret != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomDSetupAmpAndPhaseCoefficients failed\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    REAL8 PhDBMf = am * fi + bm;
    q->PhDBconst = IMRPhenomDPhase_OneFrequency(PhDBMf, pDPreComp, Rholm, Taulm) / am;

    REAL8 PhDCMf = ar * fr + br;
    q->PhDCconst = IMRPhenomDPhase_OneFrequency(PhDCMf, pDPreComp, Rholm, Taulm) / ar;

    REAL8 PhDBAMf = ai * fi + bi;
    q->PhDBAterm = IMRPhenomDPhase_OneFrequency(PhDBAMf, pDPreComp, Rholm, Taulm) / ai;
    return XLAL_SUCCESS;
}

/**
 * Define function for FD PN amplitudes
 */
COMPLEX16 IMRPhenomHMOnePointFiveSpinPN(
    REAL8 fM,
    INT4 l,
    INT4 m,
    REAL8 M1,
    REAL8 M2,
    REAL8 X1z,
    REAL8 X2z)
{

    // LLondon 2017

    // Define effective intinsic parameters
    COMPLEX16 Hlm = 0;
    REAL8 M_INPUT = M1 + M2;
    M1 = M1 / (M_INPUT);
    M2 = M2 / (M_INPUT);
    REAL8 M = M1 + M2;
    REAL8 eta = M1 * M2 / (M * M);
    REAL8 delta = sqrt(1.0 - 4 * eta);
    REAL8 Xs = 0.5 * (X1z + X2z);
    REAL8 Xa = 0.5 * (X1z - X2z);
    COMPLEX16 ans = 0;

    // Define PN parameter and realed powers
    REAL8 v = pow(M * 2.0 * LAL_PI * fM / m, 1.0 / 3.0);
    REAL8 v2 = v * v;
    REAL8 v3 = v * v2;

    // Define Leading Order Ampitude for each supported multipole
    if (l == 2 && m == 2)
    {
        // (l,m) = (2,2)
        // THIS IS LEADING ORDER
        Hlm = 1.0;
    }
    else if (l == 2 && m == 1)
    {
        // (l,m) = (2,1)
        // SPIN TERMS ADDED

        // UP TO 4PN
        REAL8 v4 = v * v3;
        Hlm = (sqrt(2.0) / 3.0) * \
            ( \
                v * delta - v2 * 1.5 * (Xa + delta * Xs) + \
                v3 * delta * ((335.0 / 672.0) + (eta * 117.0 / 56.0)
            ) \
            + \
            v4 * \
                ( \
                Xa * (3427.0 / 1344 - eta * 2101.0 / 336) + \
                delta * Xs * (3427.0 / 1344 - eta * 965 / 336) + \
                delta * (-I * 0.5 - LAL_PI - 2 * I * 0.69314718056) \
                )
            );
    }
    else if (l == 3 && m == 3)
    {
        // (l,m) = (3,3)
        // THIS IS LEADING ORDER
        Hlm = 0.75 * sqrt(5.0 / 7.0) * (v * delta);
    }
    else if (l == 3 && m == 2)
    {
        // (l,m) = (3,2)
        // NO SPIN TERMS to avoid roots
        Hlm = (1.0 / 3.0) * sqrt(5.0 / 7.0) * (v2 * (1.0 - 3.0 * eta));
    }
    else if (l == 4 && m == 4)
    {
        // (l,m) = (4,4)
        // THIS IS LEADING ORDER
        Hlm = (4.0 / 9.0) * sqrt(10.0 / 7.0) * v2 * (1.0 - 3.0 * eta);
    }
    else if (l == 4 && m == 3)
    {
        // (l,m) = (4,3)
        // NO SPIN TERMS TO ADD AT DESIRED ORDER
        Hlm = 0.75 * sqrt(3.0 / 35.0) * v3 * delta * (1.0 - 2.0 * eta);
    }
    else
    {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", l, m);
        XLAL_ERROR(XLAL_EDOM);
    }
    // Compute the final PN Amplitude at Leading Order in fM
    ans = M * M * LAL_PI * sqrt(eta * 2.0 / 3) * pow(v, -3.5) * cabs(Hlm);

    return ans;
}

/**
 * @addtogroup LALSimIMRPhenom_c
 * @{
 *
 * @name Routines for IMR Phenomenological Model "HM"
 * @{
 *
 * @author Sebastian Khan, Francesco Pannarale, Lionel London
 *
 * @brief C code for IMRPhenomHM phenomenological waveform model.
 *
 * Inspiral-merger and ringdown phenomenological, frequecny domain
 * waveform model for binary black holes systems.
 * Models not only the dominant (l,|m|) = (2,2) modes
 * but also some of the sub-domant modes too.
 * Model described in \cite PhysRevLett.120.161102.
 * The model is based on IMRPhenomD (\cite Husa:2015iqa, \cite Khan:2015jqa)
 *
 * @note The higher mode information was not calibrated to Numerical Relativity
 * simulation therefore the calibration range is inherited from PhenomD.
 *
 * @attention The model is usable outside this parameter range,
 * and in tests to date gives sensible physical results,
 * but conclusive statements on the physical fidelity of
 * the model for these parameters await comparisons against further
 * numerical-relativity simulations. For more information, see the review wiki
 * under https://git.ligo.org/waveforms/reviews/phenomhm/wikis/home
 * Also a technical document in the DCC https://dcc.ligo.org/LIGO-T1800295
 */

/**
 * Returns h+ and hx in the frequency domain.
 *
 * This function can be called in the usual sense
 * where you supply a f_min, f_max and deltaF.
 * This is the case when deltaF > 0.
 * If f_max = 0. then the default ending frequnecy is used.
 * or you can also supply a custom set of discrete
 * frequency points with which to evaluate the waveform.
 * To do this you must call this function with
 * deltaF <= 0.
 *
 */
UNUSED int XLALSimIMRPhenomHM(
    UNUSED COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
    UNUSED COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
    UNUSED REAL8Sequence *freqs,               /**< Frequency points at which to evaluate the waveform (Hz) */
    UNUSED REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    UNUSED REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    UNUSED REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    UNUSED const REAL8 distance,               /**< distance of source (m) */
    UNUSED const REAL8 inclination,            /**< inclination of source (rad) */
    UNUSED const REAL8 phiRef,                 /**< reference orbital phase (rad) */
    UNUSED const REAL8 deltaF,                 /**< Sampling frequency (Hz). To use arbitrary frequency points set deltaF <= 0. */
    UNUSED REAL8 f_ref,                        /**< Reference frequency */
    UNUSED LALDict *extraParams                /**<linked list containing the extra testing GR parameters */
)
{
    /* define and init return code for this function */
    int retcode;

    /* sanity checks on input parameters: check pointers, etc. */
    /* NOTE: a lot of checks are done in the function
     * XLALSimIMRPhenomHMGethlmModes because that can also be used
     * as a standalone function. It gets called through IMRPhenomHMCore
     * so to avoid doubling up on checks alot of the checks are done in
     * XLALSimIMRPhenomHMGethlmModes.
     */
    XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);
    XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);
    XLAL_CHECK(*hptilde == NULL, XLAL_EFAULT);
    XLAL_CHECK(*hctilde == NULL, XLAL_EFAULT);
    XLAL_CHECK(distance > 0, XLAL_EDOM, "distance must be positive.\n");

    /* main: evaluate model at given frequencies */
    retcode = 0;
    retcode = IMRPhenomHMCore(
        hptilde,
        hctilde,
        freqs,
        m1_SI,
        m2_SI,
        chi1z,
        chi2z,
        distance,
        inclination,
        phiRef,
        deltaF,
        f_ref,
        extraParams);
    XLAL_CHECK(retcode == XLAL_SUCCESS,
               XLAL_EFUNC, "IMRPhenomHMCore failed in XLALSimIMRPhenomHM.");

    /* cleanup */
    /* XLALDestroy and XLALFree any pointers. */

    return XLAL_SUCCESS;
}

/** @} */
/** @} */

/**
 * internal function that returns h+ and hx.
 * Inside this function the my bulk of the work is done
 * like the loop over frequencies.
 */
int IMRPhenomHMCore(
    UNUSED COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
    UNUSED COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
    REAL8Sequence *freqs,                      /**< GW frequecny list [Hz] */
    REAL8 m1_SI,                               /**< primary mass [kg] */
    REAL8 m2_SI,                               /**< secondary mass [kg] */
    REAL8 chi1z,                               /**< aligned spin of primary */
    REAL8 chi2z,                               /**< aligned spin of secondary */
    const REAL8 distance,                      /**< distance [m] */
    const REAL8 inclination,                   /**< inclination angle */
    const REAL8 phiRef,                        /**< orbital phase at f_ref */
    const REAL8 deltaF,                        /**< frequency spacing */
    REAL8 f_ref,                               /**< reference GW frequency */
    LALDict *extraParams                       /**< LALDict struct */
)
{
    int retcode;

    // Make a pointer to LALDict to circumvent a memory leak
    // At the end we will check if we created a LALDict in extraParams
    // and destroy it if we did.
    LALDict *extraParams_in = extraParams;

    /* evaluate all hlm modes */
    SphHarmFrequencySeries **hlms = XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms = NULL;
    retcode = 0;
    retcode = XLALSimIMRPhenomHMGethlmModes(
        hlms,
        freqs,
        m1_SI,
        m2_SI,
        chi1z,
        chi2z,
        phiRef,
        deltaF,
        f_ref,
        extraParams);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "XLALSimIMRPhenomHMGethlmModes failed");

    /* compute the frequency bounds */
    const REAL8 Mtot = (m1_SI + m2_SI) / LAL_MSUN_SI;
    PhenomHMFrequencyBoundsStorage *pHMFS;
    pHMFS = XLALMalloc(sizeof(PhenomHMFrequencyBoundsStorage));
    retcode = 0;
    retcode = init_IMRPhenomHMGet_FrequencyBounds_storage(
        pHMFS,
        freqs,
        Mtot,
        deltaF,
        f_ref);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "init_IMRPhenomHMGet_FrequencyBounds_storage failed");

    /* now we have generated all hlm modes we need to
     * multiply them with the Ylm's and sum them.
     */

    LIGOTimeGPS tC = LIGOTIMEGPSZERO; // = {0, 0}
    if (pHMFS->freq_is_uniform == 1)
    { /* 1. uniformly spaced */
        XLAL_PRINT_INFO("freq_is_uniform = True\n");
        /* coalesce at t=0 */
        /* Shift by overall length in time */
        XLAL_CHECK(
            XLALGPSAdd(&tC, -1. / deltaF),
            XLAL_EFUNC,
            "Failed to shift coalescence time to t=0,\
tried to apply shift of -1.0/deltaF with deltaF=%g.",
            deltaF);
    } /* else if 2. i.e. not uniformly spaced then we don't shift. */

    /* Allocate hptilde and hctilde */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, pHMFS->npts);
    if (!(hptilde))
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hptilde)->data->data, 0, pHMFS->npts * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);

    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, pHMFS->npts);
    if (!(hctilde))
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hctilde)->data->data, 0, pHMFS->npts * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

    /* Adding the modes to form hplus, hcross
     * - use of a function that copies XLALSimAddMode but for Fourier domain structures */
    INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */

    /* setup ModeArray */
    if (extraParams == NULL)
    {
        extraParams = XLALCreateDict();
    }
    extraParams = IMRPhenomHM_setup_mode_array(extraParams);
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams);

    /* loop over modes */
    /* at this point ModeArray should contain the list of modes
     * and therefore if NULL then something is wrong and abort.
     */
    if (ModeArray == NULL)
    {
        XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
    }
    for (UINT4 ell = 2; ell < L_MAX_PLUS_1; ell++)
    {
        for (INT4 mm = 1; mm < (INT4)ell + 1; mm++)
        { /* loop over only positive m is intentional. negative m added automatically */
            /* first check if (l,m) mode is 'activated' in the ModeArray */
            /* if activated then generate the mode, else skip this mode. */
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, mm) != 1)
            { /* skip mode */
                continue;
            } /* else: generate mode */

            COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlms, ell, mm);
            if (!(hlm))
                XLAL_ERROR(XLAL_EFUNC);

            /* We test for hypothetical m=0 modes */
            if (mm == 0)
            {
                sym = 0;
            }
            else
            {
                sym = 1;
            }
            IMRPhenomHMFDAddMode(*hptilde, *hctilde, hlm, inclination, 0., ell, mm, sym); /* The phase \Phi is set to 0 - assumes phiRef is defined as half the phase of the 22 mode h22 */
        }
    }

    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = XLALSimPhenomUtilsFDamp0(Mtot, distance);
    #pragma omp parallel for
    for (size_t i = pHMFS->ind_min; i < pHMFS->ind_max; i++)
    {
        ((*hptilde)->data->data)[i] = ((*hptilde)->data->data)[i] * amp0;
        ((*hctilde)->data->data)[i] = -1 * ((*hctilde)->data->data)[i] * amp0;
    }

    /* cleanup */
    XLALDestroyValue(ModeArray);
    LALFree(pHMFS);

  /* If extraParams was allocated in this function and not passed in
   * we need to free it to prevent a leak */
  if (extraParams && !extraParams_in) {
    XLALDestroyDict(extraParams);
  }

    return XLAL_SUCCESS;
}

/**
 * XLAL function that returns
 * a SphHarmFrequencySeries object
 * containing all the hlm modes
 * requested.
 * These have the correct relative phases between modes.
 * Note this has a similar interface to XLALSimIMRPhenomHM
 * because it is designed so it can be used independently.
 */
int XLALSimIMRPhenomHMGethlmModes(
    UNUSED SphHarmFrequencySeries **hlms, /**< [out] SphHarmFrequencySeries struct containing hlm modes */
    UNUSED REAL8Sequence *freqs,          /**< frequency sequency in Hz */
    UNUSED REAL8 m1_SI,                   /**< primary mass [kg] */
    UNUSED REAL8 m2_SI,                   /**< secondary mass [kg] */
    UNUSED REAL8 chi1z,                   /**< aligned spin of primary */
    UNUSED REAL8 chi2z,                   /**< aligned spin of secondary */
    UNUSED const REAL8 phiRef,            /**< orbital phase at f_ref */
    UNUSED const REAL8 deltaF,            /**< frequency spacing */
    UNUSED REAL8 f_ref,                   /**< reference GW frequency */
    UNUSED LALDict *extraParams           /**< LALDict struct */
)
{
    UNUSED int retcode;

    // Make a pointer to LALDict to circumvent a memory leak
    // At the end we will check if we created a LALDict in extraParams
    // and destroy it if we did.
    LALDict *extraParams_in = extraParams;

    /* sanity checks on input parameters: check pointers, etc. */

    /* Check inputs for sanity */
    XLAL_CHECK(m1_SI > 0, XLAL_EDOM, "m1 must be positive.\n");
    XLAL_CHECK(m2_SI > 0, XLAL_EDOM, "m2 must be positive.\n");
    XLAL_CHECK(fabs(chi1z) <= 1.0, XLAL_EDOM, "Aligned spin chi1z=%g \
must be <= 1 in magnitude!\n", chi1z);
    XLAL_CHECK(fabs(chi2z) <= 1.0, XLAL_EDOM, "Aligned spin chi2z=%g \
must be <= 1 in magnitude!\n", chi2z);
    XLAL_CHECK(f_ref >= 0, XLAL_EDOM, "Reference frequency must be \
positive.\n");

    /* setup ModeArray */
    if (extraParams == NULL)
        extraParams = XLALCreateDict();
    extraParams = IMRPhenomHM_setup_mode_array(extraParams);
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams);

    /* setup frequency sequency */
    REAL8Sequence *amps = NULL;
    REAL8Sequence *phases = NULL;
    REAL8Sequence *freqs_geom = NULL; /* freqs is in geometric units */

    LIGOTimeGPS tC = LIGOTIMEGPSZERO; // = {0, 0}

    /* setup PhenomHM model storage struct / structs */
    /* Compute quantities/parameters related to PhenomD only once and store them */
    PhenomHMStorage *pHM;
    pHM = XLALMalloc(sizeof(PhenomHMStorage));
    retcode = 0;
    retcode = init_PhenomHM_Storage(
        pHM,
        m1_SI,
        m2_SI,
        chi1z,
        chi2z,
        freqs,
        deltaF,
        f_ref,
        phiRef);
    XLAL_CHECK(XLAL_SUCCESS == retcode, XLAL_EFUNC, "init_PhenomHM_Storage \
failed");

    /* Two possibilities */
    if (pHM->freq_is_uniform == 1)
    { /* 1. uniformly spaced */
        XLAL_PRINT_INFO("freq_is_uniform = True\n");

        freqs = XLALCreateREAL8Sequence(pHM->npts);
        phases = XLALCreateREAL8Sequence(pHM->npts);
        amps = XLALCreateREAL8Sequence(pHM->npts);

        for (size_t i = 0; i < pHM->npts; i++)
        {                                     /* populate the frequency unitformly from zero - this is the standard
             convention we use when generating waveforms in LAL. */
            freqs->data[i] = i * pHM->deltaF; /* This is in Hz */
            phases->data[i] = 0;              /* initalise all phases to zero. */
            amps->data[i] = 0;                /* initalise all amps to zero. */
        }
        /* coalesce at t=0 */
        XLAL_CHECK(
            XLALGPSAdd(&tC, -1. / pHM->deltaF),
            XLAL_EFUNC,
            "Failed to shift coalescence time to t=0,\
tried to apply shift of -1.0/deltaF with deltaF=%g.",
            pHM->deltaF);
    }
    else if (pHM->freq_is_uniform == 0)
    { /* 2. arbitrarily space */
        XLAL_PRINT_INFO("freq_is_uniform = False\n");
        freqs = pHM->freqs; /* This is in Hz */
        phases = XLALCreateREAL8Sequence(freqs->length);
        amps = XLALCreateREAL8Sequence(freqs->length);
        for (size_t i = 0; i < pHM->npts; i++)
        {
            phases->data[i] = 0; /* initalise all phases to zero. */
            amps->data[i] = 0;   /* initalise all phases to zero. */
        }
    }
    else
    {
        XLAL_ERROR(XLAL_EDOM, "freq_is_uniform = %i and should be either 0 or 1.", pHM->freq_is_uniform);
    }

    /* PhenomD functions take geometric frequencies */
    freqs_geom = XLALCreateREAL8Sequence(pHM->npts);
    for (size_t i = 0; i < pHM->npts; i++)
    {
        freqs_geom->data[i] = XLALSimPhenomUtilsHztoMf(freqs->data[i], pHM->Mtot); /* initalise all phases to zero. */
    }

    PhenDAmpAndPhasePreComp pDPreComp22;
    retcode = IMRPhenomDSetupAmpAndPhaseCoefficients(
        &pDPreComp22,
        pHM->m1,
        pHM->m2,
        pHM->chi1z,
        pHM->chi2z,
        pHM->Rholm[2][2],
        pHM->Taulm[2][2],
        extraParams);
    if (retcode != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomDSetupAmpAndPhaseCoefficients failed\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    /* compute the reference phase shift need to align the waveform so that
     the phase is equal to phiRef at the reference frequency f_ref. */
    /* the phase shift is computed by evaluating the phase of the
    (l,m)=(2,2) mode.
    phi0 is the correction we need to add to each mode. */
    REAL8 phi_22_at_f_ref = IMRPhenomDPhase_OneFrequency(pHM->Mf_ref, pDPreComp22,
                                                         1.0, 1.0);
    REAL8 phi0 = 0.5 * phi_22_at_f_ref + phiRef;

    /* loop over modes */
    /* at this point ModeArray should contain the list of modes
     * and therefore if NULL then something is wrong and abort.
     */
    if (ModeArray == NULL)
    {
        XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
    }
    for (UINT4 ell = 2; ell < L_MAX_PLUS_1; ell++)
    {
        for (INT4 mm = 1; mm < (INT4)ell + 1; mm++)
        { /* loop over only positive m is intentional. negative m added automatically */
            /* first check if (l,m) mode is 'activated' in the ModeArray */
            /* if activated then generate the mode, else skip this mode. */
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, mm) != 1)
            { /* skip mode */
                XLAL_PRINT_INFO("SKIPPING ell = %i mm = %i\n", ell, mm);
                continue;
            } /* else: generate mode */
            XLAL_PRINT_INFO("generateing ell = %i mm = %i\n", ell, mm);

            COMPLEX16FrequencySeries *hlm = XLALCreateCOMPLEX16FrequencySeries("hlm: FD waveform", &tC, 0.0, pHM->deltaF, &lalStrainUnit, pHM->npts);
            memset(hlm->data->data, 0, pHM->npts * sizeof(COMPLEX16));
            // XLALUnitMultiply(&((*hlm)->sampleUnits), &((*hlm)->sampleUnits), &lalSecondUnit);
            retcode = 0;
            retcode = IMRPhenomHMEvaluateOnehlmMode(&hlm,
                                                    amps, phases,
                                                    freqs_geom,
                                                    pHM,
                                                    ell, mm,
                                                    phi0,
                                                    extraParams);
            XLAL_CHECK(XLAL_SUCCESS == retcode,
                       XLAL_EFUNC, "IMRPhenomHMEvaluateOnehlmMode failed");

            *hlms = XLALSphHarmFrequencySeriesAddMode(*hlms, hlm, ell, mm);

            XLALDestroyCOMPLEX16FrequencySeries(hlm);
        }
    }

    /* cleanup */
    XLALDestroyREAL8Sequence(freqs_geom);
    XLALDestroyValue(ModeArray);

    if (pHM->freq_is_uniform == 1)
    { /* 1. uniformly spaced */
        XLALDestroyREAL8Sequence(phases);
        XLALDestroyREAL8Sequence(amps);
        XLALDestroyREAL8Sequence(freqs);
    }
    else if (pHM->freq_is_uniform == 0)
    { /* 2. arbitrarily space */
        XLALDestroyREAL8Sequence(phases);
        XLALDestroyREAL8Sequence(amps);
    }
    else
    {
        XLAL_ERROR(XLAL_EDOM, "freq_is_uniform = %i and should be either 0 or 1.", pHM->freq_is_uniform);
    }

    LALFree(pHM);

  /* If extraParams was allocated in this function and not passed in
   * we need to free it to prevent a leak */
  if (extraParams && !extraParams_in) {
    XLALDestroyDict(extraParams);
  }

    return XLAL_SUCCESS;
}

/**
 * Function to compute the one hlm mode.
 * Note this is not static so that IMRPhenomPv3HM
 * can also use this function
 */
int IMRPhenomHMEvaluateOnehlmMode(
    UNUSED COMPLEX16FrequencySeries **hlm, /**< [out] One hlm mode */
    UNUSED REAL8Sequence *amps,            /**< amplitude frequency sequence */
    UNUSED REAL8Sequence *phases,          /**< phase frequency sequence */
    UNUSED REAL8Sequence *freqs_geom,      /**< dimensionless frequency sequence */
    UNUSED PhenomHMStorage *pHM,           /**< PhenomHMStorage struct */
    UNUSED UINT4 ell,                      /**< ell spherical harmonic number */
    UNUSED INT4 mm,                        /**< m spherical harmonic number */
    UNUSED REAL8 phi0,                     /**< phase shift needed to align waveform to phiRef at f_ref. */
    UNUSED LALDict *extraParams            /**< LALDict struct */
)
{
    int retcode;

    /* generate phase */
    retcode = 0;
    retcode = IMRPhenomHMPhase(
        phases,
        freqs_geom,
        pHM,
        ell, mm,
        extraParams);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "IMRPhenomHMPhase failed");

    /* generate amplitude */
    retcode = 0;
    retcode = IMRPhenomHMAmplitude(
        amps,
        freqs_geom,
        pHM,
        ell, mm,
        extraParams);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "IMRPhenomHMAmplitude failed");

    /* compute time shift */
    REAL8 t0 = IMRPhenomDComputet0(
        pHM->eta, pHM->chi1z, pHM->chi2z,
        pHM->finspin, extraParams);

    REAL8 phase_term1 = 0.0;
    REAL8 phase_term2 = 0.0;
    REAL8 Mf = 0.0; /* geometric frequency */
    /* combine together to make hlm */
    //loop over hlm COMPLEX16FrequencySeries
    for (size_t i = pHM->ind_min; i < pHM->ind_max; i++)
    {
        Mf = freqs_geom->data[i];
        phase_term1 = - t0 * (Mf - pHM->Mf_ref);
        phase_term2 = phases->data[i] - (mm * phi0);
        ((*hlm)->data->data)[i] = amps->data[i] * cexp(-I * (phase_term1 + phase_term2));
    }

    /* cleanup */

    return XLAL_SUCCESS;
}

/**
 * returns IMRPhenomHM amplitude evaluated at a set of input frequencies
 * for the l,m mode
 */
int IMRPhenomHMAmplitude(
    UNUSED REAL8Sequence *amps,       /**< [out] amplitude frequency sequence */
    UNUSED REAL8Sequence *freqs_geom, /**< dimensionless frequency sequence */
    UNUSED PhenomHMStorage *pHM,      /**< PhenomHMStorage struct */
    UNUSED UINT4 ell,                 /**< ell spherical harmonic number */
    UNUSED INT4 mm,                   /**< m spherical harmonic number */
    UNUSED LALDict *extraParams       /**< LALDict struct */
)
{
    int retcode;

    /* scale input frequencies according to PhenomHM model */
    /* LL: Map the input domain (frequencies) for this ell mm multipole
    to those appropirate for the ell=|mm| multipole */
    REAL8Sequence *freqs_amp = XLALCreateREAL8Sequence(freqs_geom->length);
    for (UINT4 i = 0; i < freqs_amp->length; i++)
    {
        freqs_amp->data[i] = IMRPhenomHMFreqDomainMap(
            freqs_geom->data[i], ell, mm, pHM, AmpFlagTrue);
    }

    /* LL: Compute the PhenomD Amplitude at the mapped l=m=2 fequencies */
    retcode = 0;
    retcode = IMRPhenomDAmpFrequencySequence(
        amps,
        freqs_amp,
        pHM->ind_min, pHM->ind_max,
        pHM->m1, pHM->m2,
        pHM->chi1z, pHM->chi2z);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
               XLAL_EFUNC, "IMRPhenomDAmpFrequencySequence failed");

    /*
    LL: Here we map the ampliude's range using two steps:
    (1) We divide by the leading order l=m=2 behavior, and then
    scale in the expected PN behavior for the multipole of interest.
    NOTE that this step is done at the mapped frequencies,
    which results in smooth behavior despite the sharp featured of the domain map.
    There are other (perhaps more intuitive) options for mapping the amplitudes,
    but these do not have the desired smooth features.
    (2) An additional scaling is needed to recover the desired PN ampitude.
    This is needed becuase only frequencies appropriate for the dominant
    quadrupole have been used thusly, so the current answer does not
    conform to PN expectations for inspiral.
    This is trikier than described here, so please give it a deeper think.
    */

    int status_in_for = XLAL_SUCCESS;
    for (UINT4 i = 0; i < freqs_amp->length; i++)
    {
        PhenomHMUsefulPowers powers_of_freq_amp;
        status_in_for = PhenomHM_init_useful_powers(
            &powers_of_freq_amp, freqs_amp->data[i]);
        if (XLAL_SUCCESS != status_in_for)
        {
            XLALPrintError("PhenomHM_init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
            retcode = status_in_for;
        }
        //new

        /* LL: Calculate the corrective factor for step #2 */

        double beta_term1 = IMRPhenomHMOnePointFiveSpinPN(
            freqs_geom->data[i],
            ell,
            mm,
            pHM->m1,
            pHM->m2,
            pHM->chi1z,
            pHM->chi2z);

        double beta=0.;
        double beta_term2=0.;
        double HMamp_term1=1.;
        double HMamp_term2=1.;
        //HACK to fix equal black hole case producing NaNs.
        //More elegant solution needed.
        if (beta_term1 == 0.){
            beta = 0.;
        } else {

            beta_term2 = IMRPhenomHMOnePointFiveSpinPN(2.0 * freqs_geom->data[i] / mm, ell, mm, pHM->m1, pHM->m2, pHM->chi1z, pHM->chi2z);
            beta = beta_term1 / beta_term2;

            /* LL: Apply steps #1 and #2 */
            HMamp_term1 = IMRPhenomHMOnePointFiveSpinPN(
                freqs_amp->data[i],
                ell,
                mm,
                pHM->m1,
                pHM->m2,
                pHM->chi1z,
                pHM->chi2z);
            HMamp_term2 = IMRPhenomHMOnePointFiveSpinPN(freqs_amp->data[i], 2, 2, pHM->m1, pHM->m2, 0.0, 0.0);

        }

        //HMamp is computed here
        amps->data[i] *= beta * HMamp_term1 / HMamp_term2;
    }

    /* cleanup */
    XLALDestroyREAL8Sequence(freqs_amp);

    return XLAL_SUCCESS;
}
/**
 * returns IMRPhenomHM phase evaluated at a set of input frequencies
 * for the l,m mode
 */
int IMRPhenomHMPhase(
    UNUSED REAL8Sequence *phases,     /**< [out] phase frequency sequence */
    UNUSED REAL8Sequence *freqs_geom, /**< dimensionless frequency sequence */
    UNUSED PhenomHMStorage *pHM,      /**< PhenomHMStorage struct */
    UNUSED UINT4 ell,                 /**< ell spherical harmonic number */
    UNUSED INT4 mm,                   /**< m spherical harmonic number */
    UNUSED LALDict *extraParams       /**< LALDict struct */
)
{
    int retcode;

    HMPhasePreComp q;
    retcode = 0;
    retcode = IMRPhenomHMPhasePreComp(&q, ell, mm, pHM, extraParams);
    if (retcode != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomHMPhasePreComp failed\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Rholm = pHM->Rholm[ell][mm];
    REAL8 Taulm = pHM->Taulm[ell][mm];

    PhenDAmpAndPhasePreComp pDPreComp;
    retcode = IMRPhenomDSetupAmpAndPhaseCoefficients(
        &pDPreComp,
        pHM->m1,
        pHM->m2,
        pHM->chi1z,
        pHM->chi2z,
        Rholm,
        Taulm,
        extraParams);
    if (retcode != XLAL_SUCCESS)
    {
        XLALPrintError("XLAL Error - IMRPhenomDSetupAmpAndPhaseCoefficients failed\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    REAL8 Mf_wf = 0.0;
    REAL8 Mf = 0.0;
    REAL8 Mfr = 0.0;
    REAL8 tmpphaseC = 0.0;
    for (UINT4 i = pHM->ind_min; i < pHM->ind_max; i++)
    {
        /* Add complex phase shift depending on 'm' mode */
        phases->data[i] = cShift[mm];
        Mf_wf = freqs_geom->data[i];
        // This if ladder is in the mathematica function HMPhase. PhenomHMDev.nb
        if (!(Mf_wf > q.fi))
        { /* in mathematica -> IMRPhenDPhaseA */
            Mf = q.ai * Mf_wf + q.bi;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pDPreComp, Rholm, Taulm) / q.ai;
        }
        else if (!(Mf_wf > q.fr))
        { /* in mathematica -> IMRPhenDPhaseB */
            Mf = q.am * Mf_wf + q.bm;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pDPreComp, Rholm, Taulm) / q.am - q.PhDBconst + q.PhDBAterm;
        }
        else if ((Mf_wf > q.fr))
        { /* in mathematica -> IMRPhenDPhaseC */
            Mfr = q.am * q.fr + q.bm;
            tmpphaseC = IMRPhenomDPhase_OneFrequency(Mfr, pDPreComp, Rholm, Taulm) / q.am - q.PhDBconst + q.PhDBAterm;
            Mf = q.ar * Mf_wf + q.br;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pDPreComp, Rholm, Taulm) / q.ar - q.PhDCconst + tmpphaseC;
        }
        else
        {
            XLALPrintError("XLAL_ERROR - should not get here - in function IMRPhenomHMPhase");
            XLAL_ERROR(XLAL_EDOM);
        }
    }

    return XLAL_SUCCESS;


}
