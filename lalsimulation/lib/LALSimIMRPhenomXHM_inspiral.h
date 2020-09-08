/*
* Copyright (C) 2019 Marta Colleoni, Cecilio García Quirós
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
//
//  LALSimIMRPhenomXHM_inspiral.h
//
//
//  Created by Marta on 06/02/2019.
//

#ifndef LALSimIMRPhenomXHM_inspiral_h
#define LALSimIMRPhenomXHM_inspiral_h

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomXHM_structs.h"
#include "LALSimIMRPhenomX_internals.h"


/********** AMPLITUDE ***************/

//Collocations Points Fits
static double IMRPhenomXHM_Insp_Amp_21_iv1(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_21_iv2(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_21_iv3(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_33_iv1(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_33_iv2(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_33_iv3(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_32_iv1(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_32_iv2(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_32_iv3(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_44_iv1(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_44_iv2(double eta, double S, double chi1, double chi2, int InspAmpFlag);
static double IMRPhenomXHM_Insp_Amp_44_iv3(double eta, double S, double chi1, double chi2, int InspAmpFlag);

/******************* Pseudo PN Amplitude coefficients ***************************/
static double IMRPhenomXHM_Inspiral_Amp_rho1(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcut3, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM);
static double IMRPhenomXHM_Inspiral_Amp_rho2(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcut3, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM);
static double IMRPhenomXHM_Inspiral_Amp_rho3(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcut3, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM);

//Ansatzs for a give frequency
static double IMRPhenomXHM_Inspiral_PNAmp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);
static double IMRPhenomXHM_Inspiral_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);
static double IMRPhenomXHM_Inspiral_Amp_NDAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);
static double IMRPhenomXHM_Inspiral_PNAmp_21Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);

//Veto amplitude
void IMRPhenomXHM_Inspiral_Amplitude_Veto(
  double *iv1, double *iv2, double *iv3,
  IMRPhenomX_UsefulPowers *powers_of_f1,
  IMRPhenomX_UsefulPowers *powers_of_f2,
  IMRPhenomX_UsefulPowers *powers_of_f3,
  IMRPhenomXHMAmpCoefficients *pAmp,
  IMRPhenomXHMWaveformStruct *pWFHM
);
static int WavyPoints(double p1, double p2, double p3);


/*********** PHASE ***************/

//Fits across parameter space
static double IMRPhenomXHM_Insp_Phase_21_lambda(double eta, double S, double chi1, double chi2, int InspPhaseFlag);
static double IMRPhenomXHM_Insp_Phase_32_lambda(double eta, double S, double chi1, double chi2, int InspPhaseFlag);
static double IMRPhenomXHM_Insp_Phase_33_lambda(double eta, double S, double chi1, double chi2, int InspPhaseFlag);
static double IMRPhenomXHM_Insp_Phase_44_lambda(double eta, double S, double chi1, double chi2, int InspPhaseFlag);
static double IMRPhenomXHM_Insp_Phase_LambdaPN(double eta, int modeInt);

//Ansatzs
static double IMRPhenomXHM_Inspiral_Phase_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMPhaseCoefficients *pPhase);
static double IMRPhenomXHM_Inspiral_Phase_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif


#endif /* LALSimIMRPhenomXHM_inspiral_h */
