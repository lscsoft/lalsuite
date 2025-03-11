#ifndef _LALSIM_IMR_PHENOMX_PNR_H
#define _LALSIM_IMR_PHENOMX_PNR_H
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

#include "LALSimIMRPhenomX_PNR_internals.h"

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

    int IMRPhenomX_PNR_GeneratePNRAngles(
        REAL8Sequence *alphaPNR,
        REAL8Sequence *betaPNR,
        REAL8Sequence *gammaPNR,
        const REAL8Sequence *freqs,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams);

    int IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies(
        REAL8Sequence *alphaPNR,
        REAL8Sequence *betaPNR,
        REAL8Sequence *gammaPNR,
        const REAL8Sequence *freqs,
        IMRPhenomXWaveformStruct *pWF_SingleSpin,
        IMRPhenomXPrecessionStruct *pPrec_SingleSpin,
        IMRPhenomX_PNR_alpha_parameters *alphaParams,
        IMRPhenomX_PNR_beta_parameters *betaParams,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams);

    int IMRPhenomX_PNR_GeneratePNRAngleInterpolants(
        IMRPhenomX_PNR_angle_spline *hm_angle_spline,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalparams);

    int IMRPhenomX_PNR_GeneratePNRAlphaForAntisymmetry(
        REAL8Sequence *alphaPNR,
        const REAL8Sequence *freqs,
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPrecessionStruct *pPrec,
        LALDict *lalParams);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMX_PNR_H */
