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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
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
    static double IMRPhenomXHM_RD_Amp_21_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_21_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_33_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_33_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_32_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_32_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_44_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_44_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_21_sigma(double eta, double S, double chi1, double chi2, int RDAmpFlag);
    static double IMRPhenomXHM_RD_Amp_33_sigma(double eta, double S, double chi1, double chi2, int RDAmpFlag);// currently constant
    static double IMRPhenomXHM_RD_Amp_32_sigma(double eta, double S, double chi1, double chi2, int RDAmpFlag);// currently constant
    static double IMRPhenomXHM_RD_Amp_44_sigma(double eta, double S, double chi1, double chi2, int RDAmpFlag);// currently constant

    //ansatz, and its derivative: analytical for no mixing and numerical for mixing
    static double IMRPhenomXHM_RD_Amp_Ansatz(double ff, IMRPhenomXHMWaveformStruct *pWF,  IMRPhenomXHMAmpCoefficients *pAmp);
    static double IMRPhenomXHM_RD_Amp_DAnsatz(double ff, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp);
    static double IMRPhenomXHM_RD_Amp_NDAnsatz(double ff, IMRPhenomXHMAmpCoefficients *pAmp,  IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXAmpCoefficients *pAmp22,  IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22);

    // Feeding the ansatz with the coefficients is how we get the final reconstruction

    //veto
    static void IMRPhenomXHM_Ringdown_Amplitude_Veto(double *pV2, double *pV3, double V4, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);


    //----PHASE-----

    // no mixing fits
    static double IMRPhenomXHM_RD_Phase_22_alpha2(double eta, double S, double chi1, double chi2, int RDPhaseFlag);
    static double IMRPhenomXHM_RD_Phase_22_alphaL(double eta, double S, double chi1, double chi2, int RDPhaseFlag);
    // 32 specific fits
    static double IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(double eta, double S, double chi1, double chi2, int RDPhaseFlag);
    static double IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(double eta, double S, double chi1, double chi2, int RDPhaseFlag);
    static double IMRPhenomXHM_Ringdown_Phase_32_p1(double eta, double S, double chi1, double chi2, int RingdownPhaseFlag);
    static double IMRPhenomXHM_Ringdown_Phase_32_p2(double eta, double S, double chi1, double chi2, int RingdownPhaseFlag);
    static double IMRPhenomXHM_Ringdown_Phase_32_p3(double eta, double S, double chi1, double chi2, int RingdownPhaseFlag);
    static double IMRPhenomXHM_Ringdown_Phase_32_p4(double eta, double S, double chi1, double chi2, int RingdownPhaseFlag);

    //ansatz
    static double IMRPhenomXHM_RD_Phase_Ansatz(double ff,IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase);
    static double IMRPhenomXHM_RD_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif


#endif 
