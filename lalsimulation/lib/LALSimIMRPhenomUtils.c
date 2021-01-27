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
 * \brief External (SWIG'd) Auxiliary functions for phenomenological model development
 *
 * Helper functions for phenomenological waveform models
 * Can be used through python SWIG wrapping
 * NOTE: The convention for naming functions in there is to use
 * the prefix 'XLALSimPhenom_'
 */

#include <lal/LALSimIMRPhenomUtils.h>
#include "LALSimIMRPhenomInternalUtils.h"
#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

// /**
//  * Example how to write an external XLAL phenom function
//  */
// void XLALSimPhenomUtilsTest(){
//     printf("Hello! I am the XLALSimPhenomUtilsTest function\n");
// }

/**
 * Convert from geometric frequency to frequency in Hz
 */
double XLALSimPhenomUtilsMftoHz(
    REAL8 Mf,       /**< Geometric frequency */
    REAL8 Mtot_Msun /**< Total mass in solar masses */
)
{
    return Mf / (LAL_MTSUN_SI * Mtot_Msun);
}

/**
 * Convert from frequency in Hz to geometric frequency
 */
double XLALSimPhenomUtilsHztoMf(
    REAL8 fHz,      /**< Frequency in Hz */
    REAL8 Mtot_Msun /**< Total mass in solar masses */
)
{
    return fHz * (LAL_MTSUN_SI * Mtot_Msun);
}

/**
 * compute the frequency domain amplitude pre-factor
 */
double XLALSimPhenomUtilsFDamp0(
    REAL8 Mtot_Msun, /**< total mass in solar masses */
    REAL8 distance   /**< distance (m) */
)
{
    return Mtot_Msun * LAL_MRSUN_SI * Mtot_Msun * LAL_MTSUN_SI / distance;
}

/**
 * Wrapper for final-spin formula based on:
 * - IMRPhenomD's FinalSpin0815() for aligned spins.
 *
 * We use their convention m1>m2
 * and put <b>all in-plane spin on the larger BH</b>.
 *
 * In the aligned limit return the FinalSpin0815 value.
 *
 * Function should reproduce
 * the function FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH
 */
REAL8 XLALSimPhenomUtilsPhenomPv2FinalSpin(
    const REAL8 m1,     /**< Mass of companion 1 (solar masses) */
    const REAL8 m2,     /**< Mass of companion 2 (solar masses) */
    const REAL8 chi1_l, /**< Aligned spin of BH 1 */
    const REAL8 chi2_l, /**< Aligned spin of BH 2 */
    const REAL8 chip   /**< Dimensionless spin in the orbital plane */
)
{
    const REAL8 M = m1 + m2;

    REAL8 af_parallel, q_factor;
    if (m1 >= m2)
    {
        q_factor = m1 / M;
        af_parallel = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1_l, chi2_l);
    }
    else
    {
        q_factor = m2 / M;
        af_parallel = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1_l, chi2_l);
    }

    REAL8 Sperp = chip * q_factor * q_factor;
    REAL8 af = copysign(1.0, af_parallel) * sqrt(Sperp * Sperp + af_parallel * af_parallel);
    return af;
}

/**
 * Helper function used in PhenomHM and PhenomPv3HM
 * Returns the final mass from the fit used in PhenomD
 */
double XLALSimPhenomUtilsIMRPhenomDFinalMass(
    REAL8 m1,    /**< mass of primary in solar masses */
    REAL8 m2,    /**< mass of secondary in solar masses */
    REAL8 chi1z, /**< aligned-spin component on primary */
    REAL8 chi2z  /**< aligned-spin component on secondary */
)
{
    int retcode = 0;
    retcode = PhenomInternal_AlignedSpinEnforcePrimaryIsm1(
        &m1,
        &m2,
        &chi1z,
        &chi2z);
    XLAL_CHECK(
        XLAL_SUCCESS == retcode,
        XLAL_EFUNC,
        "PhenomInternal_AlignedSpinEnforcePrimaryIsm1 failed");
    REAL8 Mtot = m1 + m2;
    REAL8 eta = m1 * m2 / (Mtot * Mtot);
    return (1.0 - PhenomInternal_EradRational0815(eta, chi1z, chi2z));
}

/**
 * Wrapper for final-spin formula based on:
 * - IMRPhenomD's FinalSpin0815() for aligned spins.
 *
 * We use their convention m1>m2
 * and use all in-plane spin components to determine the final spin magnitude.
 *
 * In the aligned limit return the FinalSpin0815 value.
 */
REAL8 XLALSimPhenomUtilsPhenomPv3HMFinalSpin(
    const REAL8 m1,    /**< Mass of companion 1 (solar masses) */
    const REAL8 m2,    /**< Mass of companion 2 (solar masses) */
    const REAL8 chi1x, /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi1y, /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi1z, /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2x, /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2y, /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2z  /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{

    const REAL8 M = (m1 + m2) * XLALSimPhenomUtilsIMRPhenomDFinalMass(m1, m2, chi1z, chi2z);

    REAL8 af_parallel, primary_q_factor, secondary_q_factor;
    if (m1 >= m2)
    {
        primary_q_factor = m1 / M;
        secondary_q_factor = m2 / M;
        af_parallel = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    }
    else
    {
        primary_q_factor = m2 / M;
        secondary_q_factor = m1 / M;
        af_parallel = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    }

    REAL8 S1perp = sqrt(chi1x * chi1x + chi1y * chi1y) * primary_q_factor * primary_q_factor;
    REAL8 S2perp = sqrt(chi2x * chi2x + chi2y * chi2y) * secondary_q_factor * secondary_q_factor;
    REAL8 Sperp = S1perp + S2perp;
    REAL8 af = copysign(1.0, af_parallel) * sqrt(Sperp * Sperp + af_parallel * af_parallel);
    return af;
}

/**
 * Function to compute the effective precession parameter chi_p (1408.1810)
 */
REAL8 XLALSimPhenomUtilsChiP(
    const REAL8 m1,  /**< Mass of companion 1 (solar masses) */
    const REAL8 m2,  /**< Mass of companion 2 (solar masses) */
    const REAL8 s1x, /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 s1y, /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 s2x, /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 s2y /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{
    XLAL_CHECK(m1 > 0, XLAL_EDOM, "m1 must be positive.\n");
    XLAL_CHECK(m2 > 0, XLAL_EDOM, "m2 must be positive.\n");
    XLAL_CHECK(fabs(s1x * s1x + s1y * s1y) <= 1.0, XLAL_EDOM, "|S1_perp/m1^2| must be <= 1.\n");
    XLAL_CHECK(fabs(s2x * s2x + s2y * s2y) <= 1.0, XLAL_EDOM, "|S2_perp/m2^2| must be <= 1.\n");

    const REAL8 m1_2 = m1 * m1;
    const REAL8 m2_2 = m2 * m2;

    /* Magnitude of the spin projections in the orbital plane */
    const REAL8 S1_perp = m1_2 * sqrt(s1x * s1x + s1y * s1y);
    const REAL8 S2_perp = m2_2 * sqrt(s2x * s2x + s2y * s2y);

    const REAL8 A1 = 2 + (3 * m2) / (2 * m1);
    const REAL8 A2 = 2 + (3 * m1) / (2 * m2);
    const REAL8 ASp1 = A1 * S1_perp;
    const REAL8 ASp2 = A2 * S2_perp;
    const REAL8 num = (ASp2 > ASp1) ? ASp2 : ASp1;
    const REAL8 den = (m2 > m1) ? A2 * m2_2 : A1 * m1_2;
    const REAL8 chip = num / den;

    return chip;
}

IMRPhenomPv3HMYlmStruct *XLALSimIMRPhenomPv3HMComputeYlmElements(REAL8 theta, REAL8 phi, INT4 ell)
{
    IMRPhenomPv3HMYlmStruct *p = XLALMalloc(sizeof(IMRPhenomPv3HMYlmStruct));
    memset(p, 0, sizeof(IMRPhenomPv3HMYlmStruct));

    // for each ell mode we make an array from -ell ... ell
    // of spherical harmonics and their complex conjugate
    // ylm[0] = ylm
    // ylm[1] = conj(ylm)

    int midx=0;

    if (ell == 2)
    {
        midx=0;
        for (int m = -ell ; m<=ell; m++){
            p->Ylm2[0][midx] = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, ell, m);
            p->Ylm2[1][midx] = conj(p->Ylm2[0][midx]);
            midx++;
        }

    }
    else if (ell == 3)
    {
        midx = 0;
        for (int m = -ell; m <= ell; m++)
        {
            p->Ylm3[0][midx] = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, ell, m);
            p->Ylm3[1][midx] = conj(p->Ylm3[0][midx]);
            midx++;
        }
    }
    else if (ell == 4)
    {
        midx = 0;
        for (int m = -ell; m <= ell; m++)
        {
            p->Ylm4[0][midx] = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, ell, m);
            p->Ylm4[1][midx] = conj(p->Ylm4[0][midx]);
            midx++;
        }
    }
    else
    {
        XLAL_PRINT_ERROR("ell = %i mode not possible. Returning NULL\n", ell);
        XLALFree(p);
        p = NULL;
    }

    return p;
}

// IMRPhenomPv3HMAlphaStruct *XLALSimIMRPhenomPv3HMComputeAlphaElements(UINT4 ell, REAL8 alpha)
int XLALSimIMRPhenomPv3HMComputeAlphaElements(IMRPhenomPv3HMAlphaStruct **p, UINT4 ell, REAL8 alpha)
{
    // IMRPhenomPv3HMAlphaStruct *p = XLALMalloc(sizeof(IMRPhenomPv3HMAlphaStruct));
    // memset(p, 0, sizeof(IMRPhenomPv3HMAlphaStruct));

    // for each ell mode we make an array from -ell ... ell
    // and compute the following
    // [ cexp(I * m * alpha) for m in [-ell ... ell] ]
    if (*p == NULL)
    {
        *p = XLALMalloc(sizeof(IMRPhenomPv3HMAlphaStruct));
    }


    COMPLEX16 cexp_i_alpha = cexp(+I * alpha);
    COMPLEX16 cexp_2i_alpha = cexp_i_alpha * cexp_i_alpha;
    COMPLEX16 cexp_mi_alpha = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha = cexp_mi_alpha * cexp_mi_alpha;

    COMPLEX16 cexp_3i_alpha = 0.;
    COMPLEX16 cexp_m3i_alpha = 0.;
    COMPLEX16 cexp_4i_alpha = 0.;
    COMPLEX16 cexp_m4i_alpha = 0.;

    if (ell == 2)
    {

        (*p)->alpha2[0] = cexp_m2i_alpha;
        (*p)->alpha2[1] = cexp_mi_alpha;
        (*p)->alpha2[2] = 1.0;
        (*p)->alpha2[3] = cexp_i_alpha;
        (*p)->alpha2[4] = cexp_2i_alpha;

    }
    else if (ell == 3)
    {
        cexp_3i_alpha = cexp_2i_alpha * cexp_i_alpha;
        cexp_m3i_alpha = cexp_m2i_alpha * cexp_mi_alpha;

        (*p)->alpha3[0] = cexp_m3i_alpha;
        (*p)->alpha3[1] = cexp_m2i_alpha;
        (*p)->alpha3[2] = cexp_mi_alpha;
        (*p)->alpha3[3] = 1.0;
        (*p)->alpha3[4] = cexp_i_alpha;
        (*p)->alpha3[5] = cexp_2i_alpha;
        (*p)->alpha3[6] = cexp_3i_alpha;
    }
    else if (ell == 4)
    {
        cexp_3i_alpha = cexp_2i_alpha * cexp_i_alpha;
        cexp_m3i_alpha = cexp_m2i_alpha * cexp_mi_alpha;
        cexp_4i_alpha = cexp_3i_alpha * cexp_i_alpha;
        cexp_m4i_alpha = cexp_m3i_alpha * cexp_mi_alpha;

        (*p)->alpha4[0] = cexp_m4i_alpha;
        (*p)->alpha4[1] = cexp_m3i_alpha;
        (*p)->alpha4[2] = cexp_m2i_alpha;
        (*p)->alpha4[3] = cexp_mi_alpha;
        (*p)->alpha4[4] = 1.0;
        (*p)->alpha4[5] = cexp_i_alpha;
        (*p)->alpha4[6] = cexp_2i_alpha;
        (*p)->alpha4[7] = cexp_3i_alpha;
        (*p)->alpha4[8] = cexp_4i_alpha;
    }
    else
    {
        XLAL_ERROR(XLAL_EFUNC, "ell = %i mode not possible.\n", ell);
        // XLAL_ERROR_NULL(XLAL_EINVAL, "ell = %i mode not possible.\n", ell);
    }

    return XLAL_SUCCESS;
}

/**
 * computes the wigner-d elements for -beta in batches.
 * mprime - only takes positive values as the negative values are added
 * automatically according to the symmetry
 * d^{ell}_{-m',-m} = (-1)^{m'+m} d^{ell}_{m',m}.
 */
// IMRPhenomPv3HMWignderStruct *XLALSimIMRPhenomPv3HMComputeWignerdElements(UNUSED UINT4 ell, UNUSED INT4 mprime, UNUSED REAL8 b)
int XLALSimIMRPhenomPv3HMComputeWignerdElements(IMRPhenomPv3HMWignderStruct **p, UNUSED UINT4 ell, UNUSED INT4 mprime, UNUSED REAL8 b)
{
    if (*p == NULL)
    {
        *p = XLALMalloc(sizeof(IMRPhenomPv3HMWignderStruct));
    }

    REAL8 b2 = 2. * b;
    REAL8 cosb = cos(b);
    REAL8 cos2b = cos(b2);
    REAL8 sinb = sin(b);

    REAL8 b3 = 0.;
    REAL8 cos2b_over_two = 0.;
    REAL8 cos3b = 0.;

    REAL8 cosb_over_two = cosb * ONE_OVER_TWO;
    REAL8 cos2b_fac_1 = cos2b * ONE_OVER_EIGHT + THREE_OVER_EIGHT;
    REAL8 cosb_minus_1 = cosb - 1.0;
    REAL8 cosb_plus_1 = cosb + 1.0;

    REAL8 sinb_over_two = sinb * ONE_OVER_TWO;

    if (ell == 2 && mprime == 2)
    {
        //mprime == 2
        (*p)->d22[0][0] = -cosb_over_two + cos2b_fac_1;           //m=-2
        (*p)->d22[0][1] = cosb_minus_1 * sinb_over_two;           //m=-1
        (*p)->d22[0][2] = SQRT_6 * (1. - cos2b) * ONE_OVER_EIGHT; //m=0
        (*p)->d22[0][3] = -cosb_plus_1 * sinb_over_two;           //m=1
        (*p)->d22[0][4] = cosb_over_two + cos2b_fac_1;            //m=2

        //mprime == -2
        (*p)->d22[1][0] = (*p)->d22[0][4];
        (*p)->d22[1][1] = -(*p)->d22[0][3];
        (*p)->d22[1][2] = (*p)->d22[0][2];
        (*p)->d22[1][3] = -(*p)->d22[0][1];
        (*p)->d22[1][4] = (*p)->d22[0][0];
    }
    else if (ell == 2 && mprime == 1)
    {
        REAL8 sin2b = sin(2.*b);
        cos2b_over_two = cos2b * ONE_OVER_TWO;

        //mprime == 1
        (*p)->d21[0][0] = cosb_minus_1 * sinb_over_two;           //m=-2
        (*p)->d21[0][1] = cosb_over_two - cos2b_over_two;         //m=-1
        (*p)->d21[0][2] = -SQRT_6 * sin2b * ONE_OVER_FOUR;        //m=0
        (*p)->d21[0][3] = cosb_over_two + cos2b_over_two;         //m=1
        (*p)->d21[0][4] = cosb_plus_1 * sinb_over_two;            //m=2

        //mprime == -1
        (*p)->d21[1][0] = -(*p)->d21[0][4];
        (*p)->d21[1][1] = (*p)->d21[0][3];
        (*p)->d21[1][2] = -(*p)->d21[0][2];
        (*p)->d21[1][3] = (*p)->d21[0][1];
        (*p)->d21[1][4] = -(*p)->d21[0][0];
    }
    else if (ell == 3 && mprime == 3)
    {
        b3 = 3. * b;

        cos3b = cos(b3);
        REAL8 sin2b = sin(2. * b);
        REAL8 sin3b = sin(b3);

        REAL8 FIFTEEN_OVER_32_cosb = FIFTEEN_OVER_32 * cosb;
        REAL8 THREE_OVER_16_cos2b = THREE_OVER_16 * cos2b;
        REAL8 ONE_OVER_32_cos3b = ONE_OVER_32 * cos3b;
        REAL8 THREE_OVER_16_cos2b_plus_5_OVER_16 = THREE_OVER_16_cos2b + FIVE_OVER_16;

        REAL8 FIVE_sinb = 5.0 * sinb;
        REAL8 FOUR_sin2b = 4.0 * sin2b;
        REAL8 FIFTEEN_OVER_32_cosb_plus_ONE_OVER_32_cos3b = FIFTEEN_OVER_32_cosb + ONE_OVER_32_cos3b;
        REAL8 FIVE_sinb_plus_sin3b = FIVE_sinb + sin3b;
        REAL8 FOUR_sin2b_sq = -4. * (cos2b - 1.) * ONE_OVER_TWO; // pow(sin(b), 2.0) = -(cos2b - 1.)/2.;
        REAL8 cosb_m_cos3b = cosb - cos3b;

        //mprime == 3
        (*p)->d33[0][0] = -(FIFTEEN_OVER_32_cosb_plus_ONE_OVER_32_cos3b - THREE_OVER_16_cos2b_plus_5_OVER_16); //m=-3
        (*p)->d33[0][1] = mSQRT_6_OVER_32 * (FIVE_sinb_plus_sin3b - FOUR_sin2b);                               //m=-2
        (*p)->d33[0][2] = mSQRT_15_OVER_32 * (cosb_m_cos3b - FOUR_sin2b_sq);                                   //m=-1
        (*p)->d33[0][3] = SQRT_5_OVER_16 * (-3.0 * sinb + sin3b);                                              //m=0
        (*p)->d33[0][4] = SQRT_15_OVER_32 * (cosb_m_cos3b + FOUR_sin2b_sq);                                    //m=1
        (*p)->d33[0][5] = mSQRT_6_OVER_32 * (FIVE_sinb_plus_sin3b + FOUR_sin2b);                               //m=2
        (*p)->d33[0][6] = FIFTEEN_OVER_32_cosb_plus_ONE_OVER_32_cos3b + THREE_OVER_16_cos2b_plus_5_OVER_16;    //m=3

        //mprime == -3
        (*p)->d33[1][0] = (*p)->d33[0][6];
        (*p)->d33[1][1] = -(*p)->d33[0][5];
        (*p)->d33[1][2] = (*p)->d33[0][4];
        (*p)->d33[1][3] = -(*p)->d33[0][3];
        (*p)->d33[1][4] = (*p)->d33[0][2];
        (*p)->d33[1][5] = -(*p)->d33[0][1];
        (*p)->d33[1][6] = (*p)->d33[0][0];
    }
    else if (ell == 3 && mprime == 2)
    {

        b2 = 2. * b;
        b3 = 3. * b;
        REAL8 FOUR_sin2b = 4.0 * sin(b2);
        REAL8 FIVE_sinb = 5.0 * sinb;
        REAL8 sin3b = sin(b3);
        cos3b = cos(b3);

        cos2b_over_two = cos2b * ONE_OVER_TWO;

        REAL8 FIVE_sinb_sin3b = FIVE_sinb + sin3b;

        REAL8 FIVE_OVER_16_cosb = FIVE_OVER_16 * cosb;
        REAL8 THREE_OVER_16_cos3b = THREE_OVER_16 * cos3b;
        REAL8 FIVE_OVER_16_cosb_p_THREE_OVER_16_cos3b = FIVE_OVER_16_cosb + THREE_OVER_16_cos3b;

        REAL8 THREE_sin3b = 3.0 * sin3b;
        REAL8 sinb_m_THREE_sin3b = sinb - THREE_sin3b;

        //mprime == 2
        (*p)->d32[0][0] = SQRT_6_OVER_32 * (FOUR_sin2b - FIVE_sinb_sin3b); //m=-3
        (*p)->d32[0][1] = FIVE_OVER_16_cosb_p_THREE_OVER_16_cos3b - cos2b_over_two; //m=-2
        (*p)->d32[0][2] = mSQRT_10_OVER_32 * (sinb_m_THREE_sin3b + FOUR_sin2b); //m=-1
        (*p)->d32[0][3] = SQRT_30_OVER_16 * (cosb - cos3b);            //m=0
        (*p)->d32[0][4] = SQRT_10_OVER_32 * (sinb_m_THREE_sin3b - FOUR_sin2b); //m=1
        (*p)->d32[0][5] = FIVE_OVER_16_cosb_p_THREE_OVER_16_cos3b + cos2b_over_two; //m=2
        (*p)->d32[0][6] = SQRT_6_OVER_32 * (FIVE_sinb_sin3b + FOUR_sin2b); //m=3

        //mprime == -2
        (*p)->d32[1][0] = -(*p)->d32[0][6];
        (*p)->d32[1][1] = (*p)->d32[0][5];
        (*p)->d32[1][2] = -(*p)->d32[0][4];
        (*p)->d32[1][3] = (*p)->d32[0][3];
        (*p)->d32[1][4] = -(*p)->d32[0][2];
        (*p)->d32[1][5] = (*p)->d32[0][1];
        (*p)->d32[1][6] = -(*p)->d32[0][0];
    }
    else if (ell == 4 && mprime == 4)
    {
        b2 = 2. * b;
        b3 = 3. * b;
        REAL8 b4 = 4. * b;

        cos3b = cos(b3);
        REAL8 cos4b = cos(b4);

        REAL8 sin2b = sin(b2);
        REAL8 sin3b = sin(b3);
        REAL8 sin4b = sin(b4);

        REAL8 SEVEN_OVER_16_cosb = SEVEN_OVER_16 * cosb;

        REAL8 cos3b_OVER_16 = cos3b * ONE_OVER_16;

        REAL8 sinb_14 = 14.0 * sinb;
        REAL8 sin2b_14 = 14.0 * sin2b;

        REAL8 sin3b_6  = 6.0 * sin3b;

        REAL8 sin2b_14_p_sin4b = sin2b_14 + sin4b;
        REAL8 sinb_14_p_sin3b_6 = sinb_14 + sin3b_6;

        // using sin(x)**4 = 1./8. * (3. - 4*cos(2.*x) + cos(4.*x))
        REAL8 sin_pow_4 = ONE_OVER_EIGHT * (3. - 4.*cos2b + cos4b);
        REAL8 two_sin_pow_4 = 2.0 * sin_pow_4;

        // using sin(x)**2 = -0.5 * (cos(2.*x) - 1)
        REAL8 sin_pow_2 = -ONE_OVER_TWO * (cos2b - 1.0);
        REAL8 four_sin_pow_2 = 4.0 * sin_pow_2;
        REAL8 four_sin_pow_2_m_two_sin_pow_4 = four_sin_pow_2 - two_sin_pow_4;

        // sin(x)**3 = (3*sin(x) - sin(3*x))/4.
        REAL8 four_sin_pow_3_cosb = (3.*sinb - sin3b) * cosb;

        //3.0 * sin(b) + sin(3.0 * b)
        REAL8 sin3b_m_sinb_3 = sin3b - 3. * sinb;

        //mprime == 4
        (*p)->d44[0][0] = -SEVEN_OVER_16_cosb + SEVEN_OVER_32 * cos2b - cos3b_OVER_16 + cos(b4) * ONE_OVER_128 + THIRTYFIVE_OVER_128;     //m=-4
        (*p)->d44[0][1] = mSQRT_2_OVER_64 * (sinb_14_p_sin3b_6 - sin2b_14_p_sin4b);                                                       //m=-3
        (*p)->d44[0][2] = SQRT_7_OVER_16 * (four_sin_pow_2_m_two_sin_pow_4 - cosb + cos3b);                                               //m=-2
        (*p)->d44[0][3] = SQRT_14_OVER_32 * (four_sin_pow_3_cosb + sin3b_m_sinb_3);                                                       //m=-1
        (*p)->d44[0][4] = SQRT_70_16 * sin_pow_4;                                                                                         //m=0
        (*p)->d44[0][5] = SQRT_14_OVER_32 * (sin3b_m_sinb_3 - four_sin_pow_3_cosb);                                                       //m=1
        (*p)->d44[0][6] = SQRT_7_OVER_16 * (four_sin_pow_2_m_two_sin_pow_4 + cosb - cos3b);                                               //m=2
        (*p)->d44[0][7] = mSQRT_2_OVER_64 * (sinb_14_p_sin3b_6 + sin2b_14_p_sin4b);                                                       //m=3
        (*p)->d44[0][8] = sin_pow_4 * ONE_OVER_16 - sin_pow_2 * ONE_OVER_TWO + SEVEN_OVER_16 * cosb + cos3b * ONE_OVER_16 + ONE_OVER_TWO; //m=4

        //mprime == -4
        (*p)->d44[1][0] = (*p)->d44[0][8];
        (*p)->d44[1][1] = -(*p)->d44[0][7];
        (*p)->d44[1][2] = (*p)->d44[0][6];
        (*p)->d44[1][3] = -(*p)->d44[0][5];
        (*p)->d44[1][4] = (*p)->d44[0][4];
        (*p)->d44[1][5] = -(*p)->d44[0][3];
        (*p)->d44[1][6] = (*p)->d44[0][2];
        (*p)->d44[1][7] = -(*p)->d44[0][1];
        (*p)->d44[1][8] = (*p)->d44[0][0];
    }
    else if (ell == 4 && mprime == 3)
    {
        b2 = 2. * b;
        b3 = 3. * b;
        REAL8 b4 = 4. * b;

        cos3b = cos(b3);
        REAL8 cos4b = cos(b4);

        REAL8 sin2b = sin(b2);
        REAL8 sin3b = sin(b3);
        REAL8 sin4b = sin(b4);

        REAL8 sinb_14 = 14.0 * sinb;
        REAL8 sin2b_14 = 14.0 * sin2b;

        REAL8 sin3b_6 = 6. * sin3b;
        REAL8 sin3b_3 = 3. * sin3b;

        REAL8 cosb_3 = 3. * cosb;
        REAL8 cos2b_2 = 2. * cos2b;
        REAL8 cos3b_3 = 3. * cos3b;
        REAL8 cos4b_2 = 2. * cos4b;

        REAL8 cosb_3_m_cos3b_3 = cosb_3 - cos3b_3;

        REAL8 sin2b_14_p_sin4b = sin2b_14 + sin4b;
        REAL8 sinb_14_p_sin3b_6 = sinb_14 + sin3b_6;

        REAL8 cosb_7_OVER_32 = SEVEN_OVER_32 * cosb;
        REAL8 cos2b_7_OVER_16 = SEVEN_OVER_16 * cos2b;

        REAL8 cos3b_9_OVER_32 = NINE_OVER_32 * cos3b;

        REAL8 cos4b_OVER_16 = cos4b * ONE_OVER_16;

        REAL8 cosb_7_OVER_32_p_cos3b_9_OVER_32 = cosb_7_OVER_32 + cos3b_9_OVER_32;
        REAL8 cos2b_7_OVER_16_p_cos4b_OVER_16 = cos2b_7_OVER_16 + cos4b_OVER_16;

        REAL8 sin2b_2 = 2.0 * sin2b;

        REAL8 sin2b_2_p_sin4b = sin2b_2 + sin4b;

        // sin(x)**3 = (3*sin(x) - sin(3*x))/4.
        REAL8 sin_pow_3_cosb = (3. * sinb - sin3b) * cosb * ONE_OVER_FOUR;

        //mprime == 3
        (*p)->d43[0][0] = mSQRT_2_OVER_64 * (sinb_14_p_sin3b_6 - sin2b_14_p_sin4b);           //m=-4
        (*p)->d43[0][1] = cosb_7_OVER_32_p_cos3b_9_OVER_32 - cos2b_7_OVER_16_p_cos4b_OVER_16; //m=-3
        (*p)->d43[0][2] = SQRT_14_OVER_32 * (sin3b_3 - sinb - sin2b_2_p_sin4b);               //m=-2
        (*p)->d43[0][3] = SQRT_7_OVER_32 * (cosb_3_m_cos3b_3 - cos2b_2 + cos4b_2);            //m=-1
        (*p)->d43[0][4] = -SQRT_35_OVER_4 * sin_pow_3_cosb;                                   //m=0
        (*p)->d43[0][5] = SQRT_7_OVER_32 * (cosb_3_m_cos3b_3 + cos2b_2 - cos4b_2);            //m=1
        (*p)->d43[0][6] = SQRT_14_OVER_32 * (sinb - sin3b_3 - sin2b_2_p_sin4b);               //m=2
        (*p)->d43[0][7] = cosb_7_OVER_32_p_cos3b_9_OVER_32 + cos2b_7_OVER_16_p_cos4b_OVER_16; //m=3
        (*p)->d43[0][8] = SQRT_2_OVER_64 * (sinb_14_p_sin3b_6 + sin2b_14_p_sin4b);            //m=4

        //mprime == -3
        (*p)->d43[1][0] = -(*p)->d43[0][8];
        (*p)->d43[1][1] = (*p)->d43[0][7];
        (*p)->d43[1][2] = -(*p)->d43[0][6];
        (*p)->d43[1][3] = (*p)->d43[0][5];
        (*p)->d43[1][4] = -(*p)->d43[0][4];
        (*p)->d43[1][5] = (*p)->d43[0][3];
        (*p)->d43[1][6] = -(*p)->d43[0][2];
        (*p)->d43[1][7] = (*p)->d43[0][1];
        (*p)->d43[1][8] = -(*p)->d43[0][0];
    }
    else
    {
        XLAL_ERROR(XLAL_EFUNC, "ell = %i, mprime = %i modes not possible.\n", ell, mprime);
    }

    return XLAL_SUCCESS;
};


/**
 * Hard coded expressions for relevant Wignerd matrix elements.
 * Only certain elements are coded up.
 * These were obtained using sympy and cross checked numerically
 * against the LAL XLALWignerdMatrix function
 */
int XLALSimPhenomUtilsPhenomPv3HMWignerdElement(
    REAL8 *wig_d, /**< [out] element of Wignerd Matrix */
    UINT4 ell,    /**< spherical harmonics ell mode */
    INT4 mprime,  /**< spherical harmonics m prime mode */
    INT4 mm,      /**< spherical harmonics m mode */
    REAL8 b       /**< beta angle (rad) */
)
{
    // There is alot of symmetry here so can optimise this futher
    REAL8 ans = 0.;
    REAL8 fac = 1.;

    /* to get negative mprime we use the symmetry
    d^{\ell}_{-m',-m} = (-1)^{m'+m} d^{\ell}_{m',m}
    */
    if (mprime < 0)
    {
        fac = pow(-1., mm + mprime);
        mm *= -1;
    }

    if (ell == 2)
    {
        // mprime should be the modes included in the co-precessing model
        if (abs(mprime) == 2)
        {
            if (mm == -2)
            {
                // ell =2, mprime=2, mm = -2
                ans = -cos(b) / 2.0 + cos(2.0 * b) / 8.0 + 3.0 / 8.0;
            }
            if (mm == -1)
            {
                // ell =2, mprime=2, mm = -1
                ans = (cos(b) - 1.0) * sin(b) / 2.0;
            }
            if (mm == 0)
            {
                // ell =2, mprime=2, mm = 0
                ans = sqrt(6.0) * pow(sin(b), 2.0) / 4.0;
            }
            if (mm == 1)
            {
                // ell =2, mprime=2, mm = 1
                ans = -(cos(b) + 1.0) * sin(b) / 2.0;
            }
            if (mm == 2)
            {
                // ell =2, mprime=2, mm = 2
                ans = cos(b) / 2.0 + cos(2.0 * b) / 8.0 + 3.0 / 8.0;
            }
        }
        if (abs(mprime) == 1)
        {
            if (mm == -2)
            {
                // ell =2, mprime=1, mm = -2
                ans = (cos(b) - 1.0) * sin(b) / 2.0;
            }
            if (mm == -1)
            {
                // ell =2, mprime=1, mm = -1
                ans = cos(b) / 2.0 - cos(2.0 * b) / 2.0;
            }
            if (mm == 0)
            {
                // ell =2, mprime=1, mm = 0
                ans = -sqrt(6.0) * sin(2.0 * b) / 4.0;
            }
            if (mm == 1)
            {
                // ell =2, mprime=1, mm = 1
                ans = cos(b) / 2.0 + cos(2.0 * b) / 2.0;
            }
            if (mm == 2)
            {
                // ell =2, mprime=1, mm = 2
                ans = (cos(b) + 1.0) * sin(b) / 2.0;
            }
        }
    }
    else if (ell == 3)
    {
        if (abs(mprime) == 3)
        {
            if (mm == -3)
            {
                // ell =3, mprime=3, mm = -3
                ans = -15.0 * cos(b) / 32.0 + 3.0 * cos(2.0 * b) / 16.0 - cos(3.0 * b) / 32.0 + 5.0 / 16.0;
            }
            if (mm == -2)
            {
                // ell =3, mprime=3, mm = -2
                ans = sqrt(6.0) * (-5.0 * sin(b) + 4.0 * sin(2.0 * b) - sin(3.0 * b)) / 32.0;
            }
            if (mm == -1)
            {
                // ell =3, mprime=3, mm = -1
                ans = sqrt(15.0) * (4.0 * pow(sin(b), 2.0) - cos(b) + cos(3.0 * b)) / 32.0;
            }
            if (mm == 0)
            {
                // ell =3, mprime=3, mm = 0
                ans = sqrt(5.0) * (-3.0 * sin(b) + sin(3.0 * b)) / 16.0;
            }
            if (mm == 1)
            {
                // ell =3, mprime=3, mm = 1
                ans = sqrt(15.0) * (4.0 * pow(sin(b), 2.0) + cos(b) - cos(3.0 * b)) / 32.0;
            }
            if (mm == 2)
            {
                // ell =3, mprime=3, mm = 2
                ans = -sqrt(6.0) * (5.0 * sin(b) + 4.0 * sin(2.0 * b) + sin(3.0 * b)) / 32.0;
            }
            if (mm == 3)
            {
                // ell =3, mprime=3, mm = 3
                ans = 15.0 * cos(b) / 32.0 + 3.0 * cos(2.0 * b) / 16.0 + cos(3.0 * b) / 32.0 + 5.0 / 16.0;
            }
        }
        if (abs(mprime) == 2)
        {
            if (mm == -3)
            {
                // ell =3, mprime=2, mm = -3
                ans = sqrt(6.0) * (-5.0 * sin(b) + 4.0 * sin(2.0 * b) - sin(3.0 * b)) / 32.0;
            }
            if (mm == -2)
            {
                // ell =3, mprime=2, mm = -2
                ans = 5.0 * cos(b) / 16.0 - cos(2.0 * b) / 2.0 + 3.0 * cos(3.0 * b) / 16.0;
            }
            if (mm == -1)
            {
                // ell =3, mprime=2, mm = -1
                ans = sqrt(10.0) * (-sin(b) - 4.0 * sin(2.0 * b) + 3.0 * sin(3.0 * b)) / 32.0;
            }
            if (mm == 0)
            {
                // ell =3, mprime=2, mm = 0
                ans = sqrt(30.0) * (cos(b) - cos(3.0 * b)) / 16.0;
            }
            if (mm == 1)
            {
                // ell =3, mprime=2, mm = 1
                ans = sqrt(10.0) * (sin(b) - 4.0 * sin(2.0 * b) - 3.0 * sin(3.0 * b)) / 32.0;
            }
            if (mm == 2)
            {
                // ell =3, mprime=2, mm = 2
                ans = 5.0 * cos(b) / 16.0 + cos(2.0 * b) / 2.0 + 3.0 * cos(3.0 * b) / 16.0;
            }
            if (mm == 3)
            {
                // ell =3, mprime=2, mm = 3
                ans = sqrt(6.0) * (5.0 * sin(b) + 4.0 * sin(2.0 * b) + sin(3.0 * b)) / 32.0;
            }
        }
    }
    else if (ell == 4)
    {
        if (abs(mprime) == 4)
        {
            if (mm == -4)
            {
                // ell =4, mprime=4, mm = -4
                ans = -7.0 * cos(b) / 16.0 + 7.0 * cos(2.0 * b) / 32.0 - cos(3.0 * b) / 16.0 + cos(4.0 * b) / 128.0 + 35.0 / 128.0;
            }
            if (mm == -3)
            {
                // ell =4, mprime=4, mm = -3
                ans = sqrt(2.0) * (-14.0 * sin(b) + 14.0 * sin(2.0 * b) - 6.0 * sin(3.0 * b) + sin(4 * b)) / 64.0;
            }
            if (mm == -2)
            {
                // ell =4, mprime=4, mm = -2
                ans = sqrt(7.0) * (-2.0 * pow(sin(b), 4.0) + 4.0 * pow(sin(b), 2.0) - cos(b) + cos(3.0 * b)) / 16.0;
            }
            if (mm == -1)
            {
                // ell =4, mprime=4, mm = -1
                ans = sqrt(14.0) * (4.0 * pow(sin(b), 3.0) * cos(b) - 3.0 * sin(b) + sin(3.0 * b)) / 32.0;
            }
            if (mm == 0)
            {
                // ell =4, mprime=4, mm = 0
                ans = sqrt(70.0) * pow(sin(b), 4.0) / 16.0;
            }
            if (mm == 1)
            {
                // ell =4, mprime=4, mm = 1
                ans = sqrt(14.0) * (-4.0 * pow(sin(b), 3.0) * cos(b) - 3.0 * sin(b) + sin(3.0 * b)) / 32.0;
            }
            if (mm == 2)
            {
                // ell =4, mprime=4, mm = 2
                ans = sqrt(7.0) * (-2.0 * pow(sin(b), 4.0) + 4.0 * pow(sin(b), 2.0) + cos(b) - cos(3.0 * b)) / 16.0;
            }
            if (mm == 3)
            {
                // ell =4, mprime=4, mm = 3
                ans = -sqrt(2.0) * (14.0 * sin(b) + 14.0 * sin(2.0 * b) + 6.0 * sin(3.0 * b) + sin(4.0 * b)) / 64.0;
            }
            if (mm == 4)
            {
                // ell =4, mprime=4, mm = 4
                ans = pow(sin(b), 4.0) / 16.0 - pow(sin(b), 2.0) / 2.0 + 7.0 * cos(b) / 16.0 + cos(3.0 * b) / 16.0 + 1.0 / 2.0;
            }
        }
        if (abs(mprime) == 3)
        {
            if (mm == -4)
            {
                // ell =4, mprime=3, mm = -4
                ans = sqrt(2.0) * (-14.0 * sin(b) + 14.0 * sin(2.0 * b) - 6.0 * sin(3.0 * b) + sin(4.0 * b)) / 64.0;
            }
            if (mm == -3)
            {
                // ell =4, mprime=3, mm = -3
                ans = 7.0 * cos(b) / 32.0 - 7.0 * cos(2.0 * b) / 16.0 + 9.0 * cos(3.0 * b) / 32.0 - cos(4.0 * b) / 16.0;
            }
            if (mm == -2)
            {
                // ell =4, mprime=3, mm = -2
                ans = sqrt(14.0) * (-sin(b) - 2.0 * sin(2.0 * b) + 3.0 * sin(3.0 * b) - sin(4.0 * b)) / 32.0;
            }
            if (mm == -1)
            {
                // ell =4, mprime=3, mm = -1
                ans = sqrt(7.0) * (3.0 * cos(b) - 2.0 * cos(2.0 * b) - 3.0 * cos(3.0 * b) + 2.0 * cos(4.0 * b)) / 32.0;
            }
            if (mm == 0)
            {
                // ell =4, mprime=3, mm = 0
                ans = -sqrt(35.0) * pow(sin(b), 3.0) * cos(b) / 4.0;
            }
            if (mm == 1)
            {
                // ell =4, mprime=3, mm = 1
                ans = sqrt(7.0) * (3.0 * cos(b) + 2.0 * cos(2.0 * b) - 3.0 * cos(3.0 * b) - 2.0 * cos(4.0 * b)) / 32.0;
            }
            if (mm == 2)
            {
                // ell =4, mprime=3, mm = 2
                ans = sqrt(14.0) * (sin(b) - 2.0 * sin(2.0 * b) - 3.0 * sin(3.0 * b) - sin(4.0 * b)) / 32.0;
            }
            if (mm == 3)
            {
                // ell =4, mprime=3, mm = 3
                ans = 7.0 * cos(b) / 32.0 + 7.0 * cos(2.0 * b) / 16.0 + 9.0 * cos(3.0 * b) / 32.0 + cos(4.0 * b) / 16.0;
            }
            if (mm == 4)
            {
                // ell =4, mprime=3, mm = 4
                ans = sqrt(2.0) * (14.0 * sin(b) + 14.0 * sin(2.0 * b) + 6.0 * sin(3.0 * b) + sin(4.0 * b)) / 64.0;
            }
        }
    }

    *wig_d = ans * fac;

    return XLAL_SUCCESS;
}
