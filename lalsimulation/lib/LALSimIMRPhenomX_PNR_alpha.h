#ifndef _LALSIM_IMR_PHENOMX_PNR_ALPHA_H
#define _LALSIM_IMR_PHENOMX_PNR_ALPHA_H
/*
 * Copyright (C) 2022 Cardiff University
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
 * \author Eleanor Hamilton, Sebastian Khan, Jonathan E. Thompson
 *
 */

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct tagIMRPhenomX_PNR_alpha_parameters
    {
        // alpha Ansatz terms
        REAL8 A1; /**< MR Ansatz coefficient */
        REAL8 A2; /**< MR Ansatz coefficient */
        REAL8 A3; /**< MR Ansatz coefficient */
        REAL8 A4; /**< MR Ansatz coefficient */

        // connection values
        REAL8 Mf_alpha_lower;         /**< connection frequency */
        REAL8 Mf_alpha_upper;         /**< connection frequency */
        REAL8 alpha_lower;            /**< alpha at connection frequency */
        REAL8 alpha_upper;            /**< alpha at connection frequency */
        REAL8 derivative_alpha_lower; /**< derivative at connection frequency */
        REAL8 derivative_alpha_upper; /**< derivative at connection frequency */
        REAL8 alpha_interp_0;         /**< intermediate coefficient */
        REAL8 alpha_interp_1;         /**< intermediate coefficient */
        REAL8 alpha_interp_2;         /**< intermediate coefficient */
        REAL8 alpha_interp_3;         /**< intermediate coefficient */
        REAL8 alpha_MR_offset;        /**< continuity offset between intermediate and MR */
    } IMRPhenomX_PNR_alpha_parameters;

    REAL8 IMRPhenomX_PNR_GeneratePNRAlphaAtMf(
        REAL8 Mf,
        const IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(
        REAL8 Mf,
        const IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_GetPNAlphaAtFreq(
        REAL8 Mf,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    int IMRPhenomX_PNR_precompute_alpha_coefficients(
        IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_alpha_interpolation_0(REAL8 Mf1, REAL8 Mf2, REAL8 alpha1, REAL8 alpha2, REAL8 dalpha1, REAL8 dalpha2);
    REAL8 IMRPhenomX_PNR_alpha_interpolation_1(REAL8 Mf1, REAL8 Mf2, REAL8 alpha1, REAL8 alpha2, REAL8 dalpha1, REAL8 dalpha2);
    REAL8 IMRPhenomX_PNR_alpha_interpolation_2(REAL8 Mf1, REAL8 Mf2, REAL8 alpha1, REAL8 alpha2, REAL8 dalpha1, REAL8 dalpha2);
    REAL8 IMRPhenomX_PNR_alpha_interpolation_3(REAL8 Mf1, REAL8 Mf2, REAL8 alpha1, REAL8 alpha2, REAL8 dalpha1, REAL8 dalpha2);

    REAL8 IMRPhenomX_PNR_intermediate_alpha_expression(REAL8 Mf, const IMRPhenomX_PNR_alpha_parameters *alphaParams);
    REAL8 IMRPhenomX_PNR_MR_alpha_expression(REAL8 Mf, const IMRPhenomX_PNR_alpha_parameters *alphaParams);
    int IMRPhenomX_PNR_alpha_connection_parameters(
        IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

#ifdef __cplusplus
}
#endif

#endif /*_LALSIM_IMR_PHENOMX_PNR_ALPHA_H*/
