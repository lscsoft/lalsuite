#ifndef _LALSIM_IMR_PHENOMX_PNR_BETA_H
#define _LALSIM_IMR_PHENOMX_PNR_BETA_H
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

    typedef struct tagIMRPhenomX_PNR_beta_parameters
    {
        // beta Ansatz terms
        REAL8 B0;
        REAL8 B1;
        REAL8 B2;
        REAL8 B3;
        REAL8 B4;
        REAL8 B5;

        // connection values with fixed spins
        REAL8 Mf_beta_lower;
        REAL8 Mf_beta_upper;
        REAL8 beta_lower;
        REAL8 beta_upper;
        REAL8 derivative_beta_lower;
        REAL8 derivative_beta_upper;
        REAL8 beta_rescale_1;
        REAL8 beta_rescale_2;

        // connection values with evolved spins
        REAL8 beta_interp_0;
        REAL8 beta_interp_1;
        REAL8 beta_interp_2;
        REAL8 beta_interp_3;

        // rescaled PNR beta spline
        //gsl_spline *beta_rescale_spl;
        //gsl_interp_accel *beta_rescale_acc;

    } IMRPhenomX_PNR_beta_parameters;

    REAL8 IMRPhenomX_PNR_GeneratePNRBetaAtMf(
        REAL8 Mf,
        const IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin);

    REAL8 IMRPhenomX_PNR_GenerateMergedPNRBetaAtMf(
        REAL8 Mf,
        const IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin);

    REAL8 IMRPhenomX_PNR_GeneratePNRBetaNoMR(
        REAL8 Mf,
	const IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_GenerateRingdownPNRBeta(
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_GetPNBetaAtFreq(
        REAL8 Mf,
        const IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin);

    REAL8 IMRPhenomX_PNR_GetPNBetaAtFreq_fulltwospin(
        REAL8 Mf,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
	const IMRPhenomX_PNR_beta_parameters *betaParams);

    REAL8 IMRPhenomX_PNR_PNWaveformBetaWrapper(
        REAL8 Mf,
        REAL8 MSA_beta,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_chi_calc(REAL8 m1, REAL8 L, REAL8 J0, REAL8 L0, REAL8 chi_parr, REAL8 beta);

    REAL8 IMRPhenomX_PNR_PNWaveformBeta(REAL8 Mf, REAL8 iota, REAL8 m1, REAL8 m2, REAL8 chi, REAL8 costheta);

    int IMRPhenomX_PNR_precompute_beta_coefficients(
        IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_beta_rescaling_1(REAL8 Mf, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2);
    REAL8 IMRPhenomX_PNR_beta_rescaling_2(REAL8 Mf, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2);
    REAL8 IMRPhenomX_PNR_rescale_beta_expression(REAL8 Mf, const IMRPhenomX_PNR_beta_parameters *betaParams);

    REAL8 IMRPhenomX_PNR_MR_beta_expression(REAL8 Mf, const IMRPhenomX_PNR_beta_parameters *betaParams);
    REAL8 IMRPhenomX_PNR_MR_dbeta_expression(REAL8 Mf, const IMRPhenomX_PNR_beta_parameters *betaParams);
    REAL8 IMRPhenomX_PNR_MR_ddbeta_expression(REAL8 Mf, const IMRPhenomX_PNR_beta_parameters *betaParams);
    REAL8 IMRPhenomX_PNR_MR_dddbeta_expression(REAL8 Mf, const IMRPhenomX_PNR_beta_parameters *betaParams);

    int IMRPhenomX_PNR_BetaConnectionFrequencies(
        IMRPhenomX_PNR_beta_parameters *betaParams);

    COMPLEX16 *IMRPhenomX_PNR_three_inflection_points(const IMRPhenomX_PNR_beta_parameters *betaParams);
    COMPLEX16 *IMRPhenomX_PNR_two_inflection_points(const IMRPhenomX_PNR_beta_parameters *betaParams);
    REAL8 IMRPhenomX_PNR_single_inflection_point(const IMRPhenomX_PNR_beta_parameters *betaParams);

    int IMRPhenomX_PNR_beta_connection_parameters(
        IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin);

    int IMRPhenomX_ST_PNR_beta_connection_parameters(
        IMRPhenomX_PNR_beta_parameters *betaParams,
	IMRPhenomXWaveformStruct *pWF,
	IMRPhenomXPrecessionStruct *pPrec,
	IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin);

    REAL8 IMRPhenomX_ST_PNR_beta_coeffs_0( REAL8 Mf1, REAL8 Mf2, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2 );
    REAL8 IMRPhenomX_ST_PNR_beta_coeffs_1( REAL8 Mf1, REAL8 Mf2, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2 );
    REAL8 IMRPhenomX_ST_PNR_beta_coeffs_2( REAL8 Mf1, REAL8 Mf2, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2 );
    REAL8 IMRPhenomX_ST_PNR_beta_coeffs_3( REAL8 Mf1, REAL8 Mf2, REAL8 beta1, REAL8 beta2, REAL8 dbeta1, REAL8 dbeta2 );

    REAL8 IMRPhenomX_PNR_arctan_window(REAL8 beta);

    UINT4 IMRPhenomX_PNR_AttachMRBeta(const IMRPhenomX_PNR_beta_parameters *betaParams);

#ifdef __cplusplus
}
#endif

#endif /*_LALSIM_IMR_PHENOMX_PNR_BETA_H*/
