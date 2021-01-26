#ifndef _LALSIM_IMR_PHENOMX_INSPIRAL_H
#define _LALSIM_IMR_PHENOMX_INSPIRAL_H

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

#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"


/********************************* IMRPhenomX: Amplitude Functions *********************************/
//static double IMRPhenomX_Inspiral_Amp_22_v1(double eta, double S, double dchi, double delta, int InsAmpFlag);
static double IMRPhenomX_Inspiral_Amp_22_v2(double eta, double S, double dchi, double delta, int InsAmpFlag);
static double IMRPhenomX_Inspiral_Amp_22_v3(double eta, double S, double dchi, double delta, int InsAmpFlag);
static double IMRPhenomX_Inspiral_Amp_22_v4(double eta, double S, double dchi, double delta, int InsAmpFlag);

static double IMRPhenomX_Inspiral_Amp_22_rho1(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag);
static double IMRPhenomX_Inspiral_Amp_22_rho2(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag);
static double IMRPhenomX_Inspiral_Amp_22_rho3(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag);

static double IMRPhenomX_Inspiral_Amp_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);
static double IMRPhenomX_Inspiral_Amp_22_DAnsatz(double Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);

/********************************* IMRPhenomX: Phase Functions *********************************/
static double IMRPhenomX_Inspiral_Phase_22_v3(double eta, double S, double dchi, double delta, int InspPhaseFlag);
static double IMRPhenomX_Inspiral_Phase_22_d13(double eta, double S, double dchi, double delta, int InspPhaseFlag);
static double IMRPhenomX_Inspiral_Phase_22_d23(double eta, double S, double dchi, double delta, int InspPhaseFlag);
static double IMRPhenomX_Inspiral_Phase_22_d43(double eta, double S, double dchi, double delta, int InspPhaseFlag);
static double IMRPhenomX_Inspiral_Phase_22_d53(double eta, double S, double dchi, double delta, int InspPhaseFlag);

static double IMRPhenomX_Inspiral_Phase_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase);
static double IMRPhenomX_Inspiral_Phase_22_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase);


#ifdef __cplusplus
}
#endif

#endif	// of #ifndef _LALSIM_IMR_PHENOMX_INSPIRAL_H
