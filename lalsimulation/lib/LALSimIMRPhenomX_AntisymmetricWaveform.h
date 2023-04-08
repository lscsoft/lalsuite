#ifndef _LALSIM_IMR_PHENOMX_ANTISYMMETRICWAVEFORM_H
#define _LALSIM_IMR_PHENOMX_ANTISYMMETRICWAVEFORM_H
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
 * \author Shrobana Ghosh
 *
 */

#ifdef __cplusplus
extern "C"
{
#endif

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


int IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(
    REAL8Sequence *kappa,
    const REAL8Sequence *freqs,        /**< input frequency array (Hz) */
    IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
    IMRPhenomXPrecessionStruct *pPrec /** precession struct **/
);

double GetKappa_at_frequency(REAL8 v,REAL8 delta,REAL8 Chi,REAL8 theta,REAL8 eta,double b);



int IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients(
    REAL8 *A0,
    REAL8 *phi_A0,
    REAL8 *phi_B0,
    const double MfT,
    double lina,
    double linb,
    double inveta,
    IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
    IMRPhenomXPrecessionStruct *pPrec, /** precession struct **/
    IMRPhenomXPhaseCoefficients *pPhase22 /** symmetric phase coefficients struct */
);

int IMRPhenomX_PNR_GenerateAntisymmetricWaveform(
    REAL8Sequence *antisymamp, /**< [out] Amplitude of antisymmetric (2,2) waveform */
    REAL8Sequence *antisymphase, /**< [out] Phase of antisymmetric (2,2) waveform */
    const REAL8Sequence *freqs,        /**< input frequency array (Hz) */
    IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
    IMRPhenomXPrecessionStruct *pPrec, /** precession struct **/
    LALDict *lalparams
);

#ifdef __cplusplus
}
#endif

#endif // of #ifndef _LALSIM_IMR_PHENOMX_ANTISYMMETRICWAVEFORM_H