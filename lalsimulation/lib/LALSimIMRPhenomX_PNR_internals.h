#ifndef _LALSIM_IMR_PHENOMX_PNR_INTERNALS_H
#define _LALSIM_IMR_PHENOMX_PNR_INTERNALS_H
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

#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"

#include "LALSimIMRPhenomXHM.h"
#include "LALSimIMRPhenomXHM_internals.h"

#include "LALSimIMRPhenomX_PNR_coefficients.h"
#include "LALSimIMRPhenomX_PNR_alpha.h"
#include "LALSimIMRPhenomX_PNR_beta.h"

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/LALDict.h>
#include <lal/XLALError.h>

#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>

#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

    typedef struct tagIMRPhenomX_PNR_angle_spline
    {
        gsl_spline *alpha_spline;    /**< alpha cubic spline */
        gsl_spline *beta_spline;     /**< beta cubic spline */
        gsl_spline *gamma_spline;    /**< gamma cubic spline */
        gsl_interp_accel *alpha_acc; /**< alpha cubic spline accelerator */
        gsl_interp_accel *beta_acc;  /**< beta cubic spline accelerator */
        gsl_interp_accel *gamma_acc; /**< gamma cubic spline accelerator */
    } IMRPhenomX_PNR_angle_spline;

    INT4 IMRPhenomX_PNR_GetAndSetPNRVariables(
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec);

    INT4 IMRPhenomX_PNR_PopulateStructs(
        IMRPhenomXWaveformStruct **pWF_SingleSpin,
        IMRPhenomXPrecessionStruct **pPrec_SingleSpin,
	IMRPhenomX_PNR_alpha_parameters **alphaParams,
        IMRPhenomX_PNR_beta_parameters **betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams);

    void IMRPhenomX_PNR_FreeStructs(
        IMRPhenomXWaveformStruct **pWF_SingleSpin,
        IMRPhenomXPrecessionStruct **pPrec_SingleSpin,
        IMRPhenomX_PNR_alpha_parameters **alphaParams,
        IMRPhenomX_PNR_beta_parameters **betaParams);

    INT4 IMRPhenomX_PNR_GeneratePNRGamma(
        REAL8Sequence *gamma,
        const REAL8Sequence *freqs,
        const REAL8Sequence *alpha,
        const REAL8Sequence *beta);

    INT4 IMRPhenomX_PNR_GeneratePNRGamma_FromInterpolants(
        REAL8Sequence *gamma,
        const REAL8Sequence *freqs,
        IMRPhenomX_PNR_angle_spline *ab_splines);

    REAL8 IMRPhenomX_PNR_alphadot_cosbeta(REAL8 f, IMRPhenomX_PNR_angle_spline *params);

    REAL8 IMRPhenomX_PNR_LinearFrequencyMap(REAL8 Mf, REAL8 ell, REAL8 emm, REAL8 Mf_lower, REAL8 Mf_upper, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, UINT4 INSPIRAL);
    REAL8 IMRPhenomX_PNR_LinearFrequencySlope(REAL8 emm, REAL8 Mf_lower, REAL8 Mf_upper, REAL8 Mf_RD_22, REAL8 Mf_RD_lm);
    INT4 IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies(REAL8 *Mf_low, REAL8 *Mf_high, REAL8 emmprime, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, IMRPhenomXPrecessionStruct *pPrec);
    REAL8 IMRPhenomX_PNR_HMInterpolationDeltaF(REAL8 f_min, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec);

    INT4 IMRPhenomX_PNR_CheckTwoSpin(IMRPhenomXPrecessionStruct *pPrec);

    REAL8 IMRPhenomX_PNR_AngleAtFRef(const REAL8Sequence *angle, const REAL8 f_ref, const REAL8Sequence *freqs, const REAL8 deltaF);
    REAL8 IMRPhenomX_PNR_LinearInterpolate(REAL8 a0, REAL8 a1, REAL8 f0, REAL8 f1, REAL8 feval);

    INT4 IMRPhenomX_PNR_RemapThetaJSF(
        REAL8 beta_ref,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams);

    REAL8 IMRPhenomX_PNR_GenerateEffectiveRingdownFreq(
        IMRPhenomXWaveformStruct *pWF,
        UINT4 ell,
        UINT4 emmprime,
        LALDict *lalParams /**< LAL Dictionary struct */
    );

    void IMRPhenomX_PNR_AngleParameterDebugPrint(
        IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomX_PNR_beta_parameters *betaParams);

    /* Compute window function to ensure smooth transition from tuned to PN angles outside calibration region */
    REAL8 IMRPhenomX_PNR_AnglesWindow(
        IMRPhenomXWaveformStruct *pWF,
	IMRPhenomXPrecessionStruct *pPrec
    );

    /* Compute window function which controls use of PNR coprecessing deviations. */
    REAL8 IMRPhenomX_PNR_CoprecWindow(
        IMRPhenomXWaveformStruct *pWF
    );

    /* Compute XAS phase and phase derivative a reference frequency "f_inspiral_align" */
    INT4 IMRPhenomX_PNR_SetPhaseAlignmentParams(
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec
    );

    /* Align the PNR CoPrec phase and phase derivative at
    "f_inspiral_align" by changing the effective value of
    phifRef and linb */
    void IMRPhenomX_PNR_EnforceXASPhaseAlignment(
        double* linb,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPhaseCoefficients *pPhase
    );

    /* Compute  XHM phase and phase derivative a reference frequency "f_inspiral_align" */
    INT4 IMRPhenomXHM_PNR_SetPhaseAlignmentParams(
        INT4 ell,
        INT4 emm,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams
    );

    /* Align the PNR HM CoPrec phase and phase derivative at
    "f_inspiral_align" by providing the needed phase and time shifts*/
    void IMRPhenomXHM_PNR_EnforceXHMPhaseAlignment(
        double* lina,
        double* linb,
        INT4 ell,
        INT4 emm,
        IMRPhenomXWaveformStruct *pWF,
        LALDict *lalParams
    );

    /* Function to get and or store coprec params
    into pWF and pPrec */
    INT4 IMRPhenomX_PNR_GetAndSetCoPrecParams(
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams
    );

#ifdef __cplusplus
}
#endif

#endif /*_LALSIM_IMR_PHENOMX_PNR_INTERNALS_H*/
