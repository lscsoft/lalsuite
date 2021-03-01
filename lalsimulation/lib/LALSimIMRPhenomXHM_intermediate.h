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
//
//  LALSimIMRPhenomXHM_Intermediate.h
//
//
//  Created by Marta on 06/02/2019.
//

#ifndef LALSimIMRPhenomXHM_intermediate_h
#define LALSimIMRPhenomXHM_intermediate_h


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



/*********** AMPLITUDE *****************/

//Fits int1, int2. 2 collocation points
static double IMRPhenomXHM_Inter_Amp_21_int1(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_21_int2(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_33_int1(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_33_int2(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_32_int1(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_32_int2(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_44_int1(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_44_int2(double eta, double S, double chi1, double chi2, int InterAmpFlag);

//Fits int0, dint0. Extra collocation point for EMR cases
static double IMRPhenomXHM_Inter_Amp_21_int0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_21_dint0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_33_int0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_33_dint0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_32_int0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_32_dint0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_44_int0(double eta, double S, double chi1, double chi2, int InterAmpFlag);
static double IMRPhenomXHM_Inter_Amp_44_dint0(double eta, double S, double chi1, double chi2, int InterAmpFlag);

//Coefficients of polynomial. They are feed with the some collocation points.
static double IMRPhenomXHM_Intermediate_Amp_delta0(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomXHM_Intermediate_Amp_delta1(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomXHM_Intermediate_Amp_delta2(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomXHM_Intermediate_Amp_delta3(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomXHM_Intermediate_Amp_delta4(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomXHM_Intermediate_Amp_delta5(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);

//Ansatz. Inverse of a polynomial
static double IMRPhenomXHM_Intermediate_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMAmpCoefficients *pAmp);

//Veto Functions
static void IMRPhenomXHM_Intermediate_Amplitude_Veto(double *int1, double *int2, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);
static int InsideInterval(double ftest, double fmin, double fmax);
static int CrossZeroP2(double a0, double a1, double a2, double fstart, double fend);
static int CrossZeroP3(double a0, double a1, double a2, double a3, double fstart, double fend);
static int CrossZeroP4(double a0, double a1, double a2, double a3, double a4, double fstart, double fend);
static int CrossZeroP5(double a0, double a1, double a2, double a3, double a4, double a5, double fstart, double fend);

//Functions and struct need by gsl for finding the root of a 5th order polynomial
// static double pol5 (double x, void *params);
// static double pol5_deriv(double x, void *params);
// static void pol5_fdf(double x, void *params, double *y, double *dy);
// typedef struct pol5_params
// {
//     REAL8 a0, a1, a2, a3, a4, a5;
// } pol5_params;
// static int RootPol5_finder_gsl(struct pol5_params params, double x_lo, double x_hi);

//Get the coefficients of the polynomial for a particular reconstruction indicated with IntAmpFlag
static void Update_Intermediate_Amplitude_Coefficients(IMRPhenomXHMAmpCoefficients *pAmp, int IntAmpFlag);

//Check if the polynomials cross zero and lower the order if needed in an iterative way
static void ChoosePolOrder(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp);


/************** PHASE ******************/

//Fits of the collocation points across paramter space
static double IMRPhenomXHM_Inter_Phase_21_p1(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_21_p2(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_21_p3(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_21_p4(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_21_p5(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_21_p6(double eta, double S, double chi1, double chi2, int);


static double IMRPhenomXHM_Inter_Phase_33_p1(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_33_p2(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_33_p3(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_33_p4(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_33_p5(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_33_p6(double eta, double S, double chi1, double chi2, int);

static double IMRPhenomXHM_Inter_Phase_32_p1(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_32_p2(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_32_p3(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_32_p4(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_32_p5(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_32_p6(double eta, double S, double chi1, double chi2, int);


static double IMRPhenomXHM_Inter_Phase_44_p1(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_44_p2(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_44_p3(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_44_p4(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_44_p5(double eta, double S, double chi1, double chi2, int);
static double IMRPhenomXHM_Inter_Phase_44_p6(double eta, double S, double chi1, double chi2, int);

//Ansatz
static double IMRPhenomXHM_Inter_Phase_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase);
static inline double IMRPhenomXHM_Inter_Phase_dAnsatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase);
static double IMRPhenomXHM_Inter_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif


#endif /* LALSimIMPhenomXHM_Intermediate_h */
