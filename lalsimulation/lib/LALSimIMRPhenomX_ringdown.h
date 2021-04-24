#ifndef _LALSIM_IMR_PHENOMX_RINGDOWN_H
#define _LALSIM_IMR_PHENOMX_RINGDOWN_H

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
static double IMRPhenomX_Ringdown_Amp_22_v1(double eta, double S, double dchi, double delta, int RDAmpFlag);
static double IMRPhenomX_Ringdown_Amp_22_gamma2(double eta, double S, double dchi, double delta, int RDAmpFlag);
static double IMRPhenomX_Ringdown_Amp_22_gamma3(double eta, double S, double dchi, double delta, int RDAmpFlag);

static double IMRPhenomX_Ringdown_Amp_22_PeakFrequency(double gamma2,double gamma3,double fRING,double fDAMP,int IMRPhenomXRingdownAmpVersion);

static double IMRPhenomX_Ringdown_Amp_22_Ansatz( double ff, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);
static double IMRPhenomX_Ringdown_Amp_22_DAnsatz(double ff, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);

/********************************* IMRPhenomX: Phase Functions *********************************/
static double IMRPhenomX_Ringdown_Phase_22_v4(double eta, double S, double dchi, double delta, int RDPhaseFlag);
static double IMRPhenomX_Ringdown_Phase_22_d12(double eta, double S, double dchi, double delta, int RDPhaseFlag);
static double IMRPhenomX_Ringdown_Phase_22_d24(double eta, double S, double dchi, double delta, int RDPhaseFlag);
static double IMRPhenomX_Ringdown_Phase_22_d34(double eta, double S, double dchi, double delta, int RDPhaseFlag);
static double IMRPhenomX_Ringdown_Phase_22_d54(double eta, double S, double dchi, double delta, int RDPhaseFlag);

static double IMRPhenomX_Ringdown_Phase_22_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);
static double IMRPhenomX_Ringdown_Phase_22_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);

#ifdef __cplusplus
}
#endif

#endif	// of #ifndef _LALSIM_IMR_PHENOMX_RINGDOWN_H
