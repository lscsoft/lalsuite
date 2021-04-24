#ifndef _LALSIM_IMR_PHENOMTHM_INTERNALS_H
#define _LALSIM_IMR_PHENOMTHM_INTERNALS_H

/*
 * Copyright (C) 2020 Hector Estelles
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


/*
 * \author Hector Estelles
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* LAL */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* IMRPhenomT */
#include <lal/LALSimInspiral.h>

/* Constant times employed during the code */
#define tCUT_Amp -150.0  // Ending time of inspiral region for the amplitude. Same for all modes.
#define tCUT_Freq -150.0 // Ending time of inspiral region for the frequency and phase.
#define tcpMerger -25.0 // Merger collocation point time, same for frequency and amplitude.
#define tEnd 500.0       // Maximum time for which the modes are computed. Enough to contain the non-negligible power of the ringdown radiation for all the modes.

/**** STRUCTS DEFINITION ****/

/* Waveform struct for storing binary parameters, final state parameters, and parameters needed for evaluating the calibrated fits. */
typedef struct tagIMRPhenomTWaveformStruct
{
	REAL8 m1_SI;    	/* Mass of companion 1 (kg) */
	REAL8 m2_SI;		/* Mass of companion 2 (kg) */
	REAL8 q;			/* mass ratio q>1 */
	REAL8 eta;			/* Symmetric mass ratio */
	REAL8 Mtot_SI;      /* Total mass in SI units (kg) */
	REAL8 Mtot;         /* Total mass in solar masses */
	REAL8 m1;			/* Mass 1 (larger), dimensionless (i.e. m1 \in [0,1]) */
	REAL8 m2;			/* Mass 2 (smaller), dimensionless */
	REAL8 M_sec;		/* Conversion factor, total mass in seconds */
	REAL8 delta;		/* PN symmetry coefficient */

	REAL8 fRef;			/* reference GW frequency (Hz) */
	REAL8 fmin;			/* starting GW frequency (Hz) */
	REAL8 MfRef;		/* Dimensionless reference GW frequency */
	REAL8 Mfmin;		/* Dimensionless minimum GW frequency */

	REAL8 chi1L;		/* Dimensionless aligned spin of companion 1 */
	REAL8 chi2L;		/* Dimensionless aligned spin of companion 2 */
	REAL8 Shat;			/* Dimensionless effective spin parameterisation for the callibrated fits */
	REAL8 dchi;			/* Dimensionless spin difference parameterisation for the callibrated fits */

	REAL8 Mfinal;		/* Final mass of the remnant black hole */
	REAL8 afinal;		/* Final spin of the remnant black hole */
	REAL8 afinal_prec;

	REAL8 distance;		/* Luminosity distance (m) */

	REAL8 deltaT;		/* sampling interval (s) */
	REAL8 dtM;			/* sampling interval (M) */
	REAL8 dist_sec;		/* Luminosity distance (seconds) */
	REAL8 phiRef;		/* reference orbital phase (rad) */

	REAL8 ampfac; // Amplitude conversion factor from NR to physical units
	INT4 inspVersion; // 0 (default) 22 phase/frequency reconstruction with 3 regions. Any other value, 4 regions
} IMRPhenomTWaveformStruct;

/* Struct for storing the coefficient values of the frequency and phase ansatz for the 22 mode.
   It also contains information about the minimum and reference times and the length of the waveform. */
typedef struct tagIMRPhenomTPhase22Struct
{
	REAL8 omegaPeak; // 22 frequency value at the peak amplitude time of the 22 mode (t=0)

	/* Ringdown ansatz phenomenological coefficients */
	REAL8 c1;
	REAL8 c1_prec;
	REAL8 c2;
	REAL8 c3;
	REAL8 c4;

	/* PN coefficients of TaylorT3 */
	REAL8 omega1PN;
	REAL8 omega1halfPN;
	REAL8 omega2PN;
	REAL8 omega2halfPN;
	REAL8 omega3PN;
	REAL8 omega3halfPN;

	/* Phenomenological coefficients of inspiral ansatz */
	REAL8 omegaInspC1;
	REAL8 omegaInspC2;
	REAL8 omegaInspC3;
	REAL8 omegaInspC4;
	REAL8 omegaInspC5;
	REAL8 omegaInspC6;

	/* Phenomenological coefficients of merger ansatz */
	REAL8 omegaMergerC1;
	REAL8 omegaMergerC2;
	REAL8 omegaMergerC3;

	REAL8 alpha1RD; // Angular damping frequency of the nlm = 122 QNM
	REAL8 alpha1RD_prec;
	REAL8 domegaPeak;	// Frequency derivative at the peak amplitude time
	REAL8 omegaRING;	// Angular ringdown frequency
	REAL8 omegaRING_prec;

	REAL8 MfRef;  // Dimensionless reference frequency
	REAL8 Mfmin;  // Dimensionless minimum frequency	

	REAL8 tRef;		// Reference time
	REAL8 tmin;		// Minimum time
	REAL8 tminSec;	// Minimum time in seconds
	size_t wflength;	// Length of the waveform (in number of points)
	size_t wflength_insp_early;	// Length of the early inspiral waveform (in number of points) (needed if 4 regions are requested)
	size_t wflength_insp_late;  // Length of the late inspiral waveform (in number of points) (needed if 4 regions are requested)

	REAL8 phOffInsp;	// phase offset of inspiral phase (needed if 4 regions are requested)
	REAL8 phOffMerger; // phase offset of merger phase
	REAL8 phOffRD;     // phase offset of ringdown phase

	REAL8 tCut22;  // Inspiral-merger boundary time for the 22 mode phase/frequency
	REAL8 tEarly;  // Early inspiral-late inspiral boundary time for the 22 mode phase/frequency (needed if 4 regions are requested)
	REAL8 tt0;     // Calibrated t0 parameter of TaylorT3, needed for theta=0.33 collocation point and for early inspiral region if requested.

	REAL8 dtM;     // Dimensionless time step
	REAL8 EulerRDslope; // Slope of the analytical ringdown approximation for the precessing alpha angle. FIXME: check if its needed.

} IMRPhenomTPhase22Struct;

/* Struct for storing the coefficient values of the frequency and phase merger-ringdown ansatz for the lm mode. */
typedef struct tagIMRPhenomTHMPhaseStruct
{
	UINT4 emm;  // mode number
	REAL8 omegaPeak; // lm frequency value at the peak amplitude time of the lm mode

	/* Ringdown ansatz phenomenological coefficients */
	REAL8 c1;
	REAL8 c1_prec;
	REAL8 c2;
	REAL8 c3;
	REAL8 c4;

	/* Phenomenological coefficients of merger ansatz */
	REAL8 omegaMergerC1;
	REAL8 omegaMergerC2;
	REAL8 omegaMergerC3;

	REAL8 alpha1RD;	// Angular damping frequency of the nlm = 1lm QNM
	REAL8 alpha1RD_prec;
	REAL8 alpha2RD; // Angular damping frequency of the nlm = 2lm QNM
	REAL8 alpha21RD; // Damping frequency difference between n=2 and n=1 QNM
	
	REAL8 domegaPeak;	// Frequency derivative at the peak amplitude time
	REAL8 omegaRING;	// Angular ringdown frequency
	REAL8 omegaRING_prec;

	REAL8 phOffMerger;	// phase offset of merger phase
	REAL8 phOffRD;		// phase offset of ringdown phase
	
} IMRPhenomTHMPhaseStruct;

/* Struct for storing the coefficient values of the amplitude ansatz for the lm mode. */
typedef struct tahIMRPhenomTHMAmpStruct
{
	/* PN amplitude coefficients */
	REAL8 fac0;
	REAL8 ampN;
	REAL8 amp0halfPNreal;
	REAL8 amp0halfPNimag;
	REAL8 amp1PNreal;
	REAL8 amp1PNimag;
	REAL8 amp1halfPNreal;
	REAL8 amp1halfPNimag;
	REAL8 amp2PNreal;
	REAL8 amp2PNimag;
	REAL8 amp2halfPNreal;
	REAL8 amp2halfPNimag;
	REAL8 amp3PNreal;
	REAL8 amp3PNimag;
	REAL8 amp3halfPNreal;
	REAL8 amp3halfPNimag;
	REAL8 amplog;

	REAL8 tshift; // Peak amplitude time of the lm mode

	/* Phenomenological coefficients of inspiral ansatz */
	REAL8 inspC1;
	REAL8 inspC2;
	REAL8 inspC3;

	/* Phenomenological coefficients of merger ansatz */
	REAL8 mergerC1;
	REAL8 mergerC2;
	REAL8 mergerC3;
	REAL8 mergerC4;

	/* Damping frequencies */
	REAL8 alpha1RD;
	REAL8 alpha2RD;
	REAL8 alpha21RD;

	REAL8 alpha1RD_prec;
	REAL8 alpha2RD_prec;
	REAL8 alpha21RD_prec;

	/* Phenomenological coefficients of ringdown ansatz */
	REAL8 c1;
	REAL8 c2;
	REAL8 c3;
	REAL8 c4;

	REAL8 c1_prec;
	REAL8 c2_prec;
	REAL8 c4_prec;

	/* Complex inspiral amplitude phase/frequenct contribution at the inspiral-merger boundary */
	REAL8 omegaCutPNAMP;
	REAL8 phiCutPNAMP;
} IMRPhenomTHMAmpStruct;

/* Wrapper struct for GSL Root Finder */
struct FindRootParams{
	double f0;
	IMRPhenomTWaveformStruct *wf;
	IMRPhenomTPhase22Struct *pPhase;
};

/**** FUNCTIONS FOR POPULATING THE STRUCTS *****/

double GetTimeOfFreq(double t, void *params); // Wrapper function for GSL Root Finder

int IMRPhenomTSetWaveformVariables(
	IMRPhenomTWaveformStruct *wf,
	const REAL8 m1_SI,
  	const REAL8 m2_SI,
  	const REAL8 chi1L_In,
  	const REAL8 chi2L_In,
  	const REAL8 distance,
  	const REAL8 deltaT,
  	const REAL8 fmin,
  	const REAL8 fRef,
  	const REAL8 phiRef,
  	LALDict *lalParams
);

/* Functions to populate phase and amplitude structs */
int IMRPhenomTSetPhase22Coefficients(IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf);
int IMRPhenomTSetHMPhaseCoefficients(int l, int m, IMRPhenomTHMPhaseStruct *pPhaseHM, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTHMAmpStruct *pAmp, IMRPhenomTWaveformStruct *wf);
int IMRPhenomTSetHMAmplitudeCoefficients(int l, int m, IMRPhenomTHMAmpStruct *pAmp, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf);

/********************************* 22 Frequency ansatzs *********************************/

double IMRPhenomTTaylorT3(REAL8 theta, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiralOmegaAnsatz22(REAL8 theta, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTMergerOmegaAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTRDOmegaAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase);

/********************************* 22 Phase ansatzs *********************************/

double IMRPhenomTInspiralPhaseTaylorT3(REAL8 thetabar, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiralPhaseAnsatz22(REAL8 t, REAL8 theta, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTMergerPhaseAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTRDPhaseAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase);

/********************************* HM Frequency ansatzs *********************************/

double IMRPhenomTMergerOmegaAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase);
double IMRPhenomTRDOmegaAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase);

/********************************* HM Phase ansatzs *********************************/

double IMRPhenomTMergerPhaseAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase);
double IMRPhenomTRDPhaseAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase);

/********************************* HM Amplitude ansatzs *********************************/

COMPLEX16 IMRPhenomTInspiralAmpAnsatzHM(REAL8 t, IMRPhenomTHMAmpStruct *pAmp);
double IMRPhenomTMergerAmpAnsatzHM(REAL8 t, IMRPhenomTHMAmpStruct *pAmp);
double IMRPhenomTRDAmpAnsatzHM(REAL8 t, IMRPhenomTHMAmpStruct *pAmp);

/* Wrapper for complex inspiral amplitude phase */
double ComplexAmpOrientation(REAL8 xref, IMRPhenomTHMAmpStruct *pAmp);

/***************************** Piecewise functions for 22 phase and frequency ********/

double IMRPhenomTomega22(REAL8 t, REAL8 theta, IMRPhenomTWaveformStruct *pWF, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTPhase22(REAL8 t, REAL8 thetabar, IMRPhenomTWaveformStruct *pWF, IMRPhenomTPhase22Struct *pPhase);

/***************************** Piecewise functions for lm phase and amplitude ********/

double IMRPhenomTHMPhase(
  	REAL8 t,
  	REAL8 phiInsp,
  	IMRPhenomTHMPhaseStruct *pPhaseHM,
  	IMRPhenomTHMAmpStruct *pAmpHM
);

COMPLEX16 IMRPhenomTHMAmp(REAL8 t, REAL8 w, IMRPhenomTHMAmpStruct *pAmpHM);

REAL8 GetEulerSlope(REAL8 af, REAL8 mf);

#ifdef __cplusplus
}
#endif

#endif
