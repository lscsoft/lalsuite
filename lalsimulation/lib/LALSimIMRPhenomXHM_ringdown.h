/*
 * Copyright (C) 2019 Marta Colleoni, Cecilio Garcia Quiros
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
 *
 */


#ifndef LALSIMIMRPHENOMXHM_RINGDOWN_H
#define LALSIMIMRPHENOMXHM_RINGDOWN_H

#ifdef __cplusplus
extern "C" {
#endif

    #ifdef __GNUC__
    #define UNUSED __attribute__ ((unused))
    #else
    #define UNUSED
    #endif

#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"

#include "LALSimIMRPhenomXHM_structs.h"


//----AMPLITUDE----

//Fits of the ringdown coefficients over parameter space
static double IMRPhenomXHM_RD_Amp_21_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_21_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_44_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_44_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_21_sigma(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_sigma(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);// currently constant
static double IMRPhenomXHM_RD_Amp_32_sigma(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);// currently constant
static double IMRPhenomXHM_RD_Amp_44_sigma(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);// currently constant
static double IMRPhenomXHM_RD_Amp_21_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_21_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_21_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_33_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_44_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_44_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_44_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_rdaux1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);
static double IMRPhenomXHM_RD_Amp_32_rdaux2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag);

/* End of Amp Parameter Space Fits */

// Ringdown coefficients from collocation points
static void IMRPhenomXHM_RD_Amp_Coefficients(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);
static void IMRPhenomXHM_RDAux_Amp_Coefficients(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);

//ansatz, and its derivative: analytical for no mixing and numerical for mixing
static double IMRPhenomXHM_RD_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWF,  IMRPhenomXHMAmpCoefficients *pAmp);
static double IMRPhenomXHM_RD_Amp_DAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp);
static double IMRPhenomXHM_RD_Amp_NDAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp,  IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXAmpCoefficients *pAmp22,  IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22);

// Feeding the ansatz with the coefficients is how we get the final reconstruction

//veto
static void IMRPhenomXHM_Ringdown_Amplitude_Veto(double *pV2, double *pV3, double V4, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);


//----PHASE-----

// no mixing fits
static double IMRPhenomXHM_RD_Phase_22_alpha2(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);
static double IMRPhenomXHM_RD_Phase_22_alphaL(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);
// 32 specific fits
static double IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p1(IMRPhenomXWaveformStruct *pWF, int RingdownPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p2(IMRPhenomXWaveformStruct *pWF, int RingdownPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p3(IMRPhenomXWaveformStruct *pWF, int RingdownPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p4(IMRPhenomXWaveformStruct *pWF, int RingdownPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p5(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);
static double IMRPhenomXHM_RD_Phase_32_p5(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag);

/* End of Phase Parameter Space Fits */

//ansatz
static double IMRPhenomXHM_RD_Phase_Ansatz(double ff,IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase);
static double IMRPhenomXHM_RD_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif


#endif
