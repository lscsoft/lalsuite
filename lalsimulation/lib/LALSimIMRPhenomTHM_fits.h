#ifndef _LALSIM_IMR_PHENOMT_FITS_H
#define _LALSIM_IMR_PHENOMT_FITS_H

/*
 * Copyright (C) 2020 Hector Estelles
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
 * \author Hector Estelles
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRPhenomXUtilities.h"


/********************************* IMRPhenomT 22 Frequency Fits *********************************/

/*static double IMRPhenomT_MECOTime(double eta, double S, double dchi, double delta);/*FIXME: Promote some interesting quantities to XLAL functions*/

static double IMRPhenomT_Inspiral_TaylorT3_t0(double eta, double S, double dchi, double delta); // theta = 0.45

static double IMRPhenomT_Inspiral_Freq_CP1_22(double eta, double S, double dchi, double delta); // theta = 0.45
static double IMRPhenomT_Inspiral_Freq_CP2_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Freq_CP3_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Freq_CP4_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Freq_CP5_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Freq_CP1_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakFrequency_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D2_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D3_22(double eta, double S, double dchi, double delta);

/********************************* IMRPhenomT 22 Amplitude Fits *********************************/

static double IMRPhenomT_Inspiral_Amp_CP1_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP2_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP3_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Amp_CP1_22(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakAmp_22(double eta, double S, double dchi);
static double IMRPhenomT_RD_Amp_C3_22(double eta, double S);

/********************************* IMRPhenomT 22 QNM Fits *********************************/

static double evaluate_QNMfit_fring22(double finalDimlessSpin);

static double evaluate_QNMfit_fdamp22(double finalDimlessSpin);

static double evaluate_QNMfit_fdamp22n2(double finalDimlessSpin);

/********************************* IMRPhenomTHM Fits *********************************/

static double IMRPhenomT_Merger_Freq_CP1_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakFrequency_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP1_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP2_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP3_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Amp_CP1_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakAmp_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Amp_C3_21(double eta, double S, double dchi);
static double IMRPhenomT_RD_Freq_D2_21(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D3_21(double eta, double S, double dchi, double delta);

static double IMRPhenomT_Merger_Freq_CP1_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakFrequency_33(double eta, double S, double dchi);
static double IMRPhenomT_Inspiral_Amp_CP1_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP2_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP3_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Amp_CP1_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakAmp_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Amp_C3_33(double eta, double S);
static double IMRPhenomT_RD_Freq_D2_33(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D3_33(double eta, double S, double dchi, double delta);

static double IMRPhenomT_Merger_Freq_CP1_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakFrequency_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP1_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP2_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP3_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Amp_CP1_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakAmp_44(double eta, double S, double dchi), double delta;
static double IMRPhenomT_RD_Amp_C3_44(double eta, double S);
static double IMRPhenomT_RD_Freq_D2_44(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D3_44(double eta, double S, double dchi, double delta);

static double IMRPhenomT_Merger_Freq_CP1_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakFrequency_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP1_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP2_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Inspiral_Amp_CP3_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_Merger_Amp_CP1_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_PeakAmp_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Amp_C3_55(double eta, double S, double dchi);
static double IMRPhenomT_RD_Freq_D2_55(double eta, double S, double dchi, double delta);
static double IMRPhenomT_RD_Freq_D3_55(double eta, double S, double dchi, double delta);

/********************************* IMRPhenomTHM QNM Fits *********************************/

static double evaluate_QNMfit_fring21(double finalDimlessSpin);
static double evaluate_QNMfit_fring33(double finalDimlessSpin);
static double evaluate_QNMfit_fring44(double finalDimlessSpin);
static double evaluate_QNMfit_fring55(double finalDimlessSpin);

static double evaluate_QNMfit_fdamp21(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp33(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp44(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp55(double finalDimlessSpin);

static double evaluate_QNMfit_fdamp21n2(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp33n2(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp44n2(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp55n2(double finalDimlessSpin);

/***************************** IMRPhenomTHM Time Shifts ***********************/

static double IMRPhenomT_tshift_21(double eta, double S, double dchi);
static double IMRPhenomT_tshift_33(double eta, double S);
static double IMRPhenomT_tshift_44(double eta, double S);
static double IMRPhenomT_tshift_55(double eta, double S);

#ifdef __cplusplus
}
#endif

#endif	// of #ifndef _LALSIM_IMR_PHENOMT_INSPIRAL_H
