/*
 * Copyright (C) 2011 J. Creighton
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
 */

/**
 * @defgroup LALSimNoise_h Header LALSimNoise.h
 * @ingroup lalsimulation_noise
 * @author Jolien Creighton
 * @brief Routines for simulating gravitational wave detector noise.
 *
 * @{
 * @defgroup LALSimNoise_c    Module LALSimNoise.c
 * @defgroup LALSimNoisePSD_c Module LALSimNoisePSD.c
 * @}
 */

#ifndef _LALSIMNOISE_H
#define _LALSIMNOISE_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <gsl/gsl_rng.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#ifdef SWIG // SWIG interface directives
SWIGLAL(FUNCTION_POINTER(XLALSimNoisePSDiLIGOSRD, XLALSimNoisePSDiLIGOModel,
	XLALSimNoisePSDeLIGOModel, XLALSimNoisePSDVirgo, XLALSimNoisePSDGEO,
	XLALSimNoisePSDTAMA, XLALSimNoisePSDaLIGONoSRMLowPower,
	XLALSimNoisePSDaLIGONoSRMHighPower,
	XLALSimNoisePSDaLIGOZeroDetLowPower,
	XLALSimNoisePSDaLIGOZeroDetHighPower, XLALSimNoisePSDaLIGONSNSOpt,
	XLALSimNoisePSDaLIGOBHBH20Deg, XLALSimNoisePSDaLIGOHighFrequency,
	XLALSimNoisePSDKAGRA, XLALSimNoisePSDAdvVirgo));
#endif


/*
 * NOISE GENERATION ROUTINES
 * in module LALSimNoise.c
 */


int XLALSimNoise(REAL8TimeSeries *s, size_t stride, REAL8FrequencySeries *psd, gsl_rng *rng);


/*
 * PSD GENERATION FUNCTIONS
 * in module LALSimNoisePSD.c
 */


/*
 * FUNCTIONS TO GENERATE COMPONENT NOISE PSD
 */

double XLALSimNoisePSDSeismic(double f, double L, double f_pend, double f_stack, double n_stack);
double XLALSimNoisePSDSuspTherm(double f, double L, double M, double T, double f0, double Q);
double XLALSimNoisePSDMirrorTherm(double f, double L, double M, double T, double f0, double Q);
double XLALSimNoisePSDShot(double f, double P_BS, double lambda, double L, double finesse, double eta);
double XLALSimNoisePSDQuantum(double f, double I0, double lambda, double L, double M, double A, double A_BS, double T_ITM, double T_PRM, double T_SRM, double ds, double zeta, double eta);

/*
 * NOISE PSD ROUTINES FOR FIRST GENERATION DETECTORS
 */

double XLALSimNoisePSDiLIGOSRD(double f);
double XLALSimNoisePSDiLIGOSeismic(double f);
double XLALSimNoisePSDiLIGOThermal(double f);
double XLALSimNoisePSDiLIGOShot(double f);
double XLALSimNoisePSDeLIGOShot(double f);

double XLALSimNoisePSDiLIGOModel(double f);
double XLALSimNoisePSDeLIGOModel(double f);
double XLALSimNoisePSDVirgo(double f);
double XLALSimNoisePSDGEO(double f);
double XLALSimNoisePSDGEOHF(double f);
double XLALSimNoisePSDTAMA(double f);

/*
 * NOISE PSD ROUTINES FOR SECOND GENERATION DETECTORS
 */

double XLALSimNoisePSDaLIGOThermal(double f);
double XLALSimNoisePSDaLIGOQuantumNoSRMLowPower(double f);
double XLALSimNoisePSDaLIGOQuantumNoSRMHighPower(double f);
double XLALSimNoisePSDaLIGOQuantumZeroDetLowPower(double f);
double XLALSimNoisePSDaLIGOQuantumZeroDetHighPower(double f);
double XLALSimNoisePSDaLIGOQuantumNSNSOpt(double f);
double XLALSimNoisePSDaLIGOQuantumBHBH20Deg(double f);
double XLALSimNoisePSDaLIGOQuantumHighFrequency(double f);

double XLALSimNoisePSDaLIGONoSRMLowPower(double f);
double XLALSimNoisePSDaLIGONoSRMHighPower(double f);
double XLALSimNoisePSDaLIGOZeroDetLowPower(double f);
double XLALSimNoisePSDaLIGOZeroDetHighPower(double f);
double XLALSimNoisePSDaLIGONSNSOpt(double f);
double XLALSimNoisePSDaLIGOBHBH20Deg(double f);
double XLALSimNoisePSDaLIGOHighFrequency(double f);
double XLALSimNoisePSDKAGRA(double f);
double XLALSimNoisePSDAdvVirgo(double f);

/*
 * NOISE PSD UTILITY ROUTINES
 */

int XLALSimNoisePSD(REAL8FrequencySeries *psd, double flow, double (*psdfunc)(double));
int XLALSimNoisePSDFromFile(REAL8FrequencySeries *psd, double flow, const char *fname);

/*
 * NOISE PSDs FROM LIGO-T0900288
 */

int XLALSimNoisePSDaLIGONoSRMLowPowerGWINC(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOZeroDetLowPowerGWINC(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOZeroDetHighPowerGWINC(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGONSNSOptGWINC(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOBHBH20DegGWINC(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOHighFrequencyGWINC(REAL8FrequencySeries *psd, double flow);

/*
 * NOISE PSDs FROM LIGO-P1200087
 */

int XLALSimNoisePSDaLIGOEarlyLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOEarlyHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOMidLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOMidHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOLateLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOLateHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGODesignSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDaLIGOBNSOptimizedSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVEarlyLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVEarlyHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVMidLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVMidHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVLateLowSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVLateHighSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVDesignSensitivityP1200087(REAL8FrequencySeries *psd, double flow);
int XLALSimNoisePSDAdVBNSOptimizedSensitivityP1200087(REAL8FrequencySeries *psd, double flow);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMNOISE_H */
