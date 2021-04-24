#ifndef _LALSIM_IMR_PHENOMX_INTERMEDIATE_H
#define _LALSIM_IMR_PHENOMX_INTERMEDIATE_H

/*
 * Copyright (C) 2018 Geraint Pratten
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
 * \author Geraint Pratten
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"


/********************************* IMRPhenomX: Amplitude Functions *********************************/
static double IMRPhenomX_Intermediate_Amp_22_vA(double eta, double S, double dchi, double delta, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_v2(double eta, double S, double dchi, double delta, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_v3(double eta, double S, double dchi, double delta, int IntAmpFlag);

static double IMRPhenomX_Intermediate_Amp_22_delta0(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_delta1(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_delta2(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_delta3(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_delta4(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);
static double IMRPhenomX_Intermediate_Amp_22_delta5(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag);

static double IMRPhenomX_Intermediate_Amp_22_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);

/********************************* IMRPhenomX: Phase Functions *********************************/
static double IMRPhenomX_Intermediate_Phase_22_v2( double eta, double S, double dchi, double delta, int IntPhaseFlag);
UNUSED static double IMRPhenomX_Intermediate_Phase_22_v3( double eta, double S, double dchi, double delta, int IntPhaseFlag);
static double IMRPhenomX_Intermediate_Phase_22_v2mRDv4(double eta, double S, double dchi, double delta, int IntPhaseFlag);
static double IMRPhenomX_Intermediate_Phase_22_v3mRDv4(double eta, double S, double dchi, double delta, int IntPhaseFlag);
static double IMRPhenomX_Intermediate_Phase_22_d43(double eta, double S, double dchi, double delta, int IntPhaseFlag);

static double IMRPhenomX_Intermediate_Phase_22_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);
static double IMRPhenomX_Intermediate_Phase_22_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif

#endif	// of #ifndef _LALSIM_IMR_PHENOMX_INTERMEDIATE_H
