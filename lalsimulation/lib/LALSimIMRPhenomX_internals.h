#ifndef _LALSIM_IMR_PHENOMX_INTERNALS_H
#define _LALSIM_IMR_PHENOMX_INTERNALS_H

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


/*
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

/* Inherited from IMRPhenomD */
#define N_MAX_COLLOCATION_POINTS_PHASE_RD 5
#define N_MAX_COLLOCATION_POINTS_PHASE_INT 5
#define N_MAX_COLLOCATION_POINTS_PHASE_INS 6

#define N_MAX_COLLOCATION_POINTS_AMP_RD 5
#define N_MAX_COLLOCATION_POINTS_AMP_INT 5
#define N_MAX_COLLOCATION_POINTS_AMP_INS 5

/* Standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* LAL */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

/* IMRPhenomX */
#include <lal/LALSimInspiral.h>

/* ********************** CACHED VARIABLES ********************* */
/*
	Cached, recurring coefficients that occur in IMRPhenomX.
	Must be frequency independent.
*/
typedef struct tagIMRPhenomXWaveformStruct
{
	/* Debug flag */
	INT4 debug;

	/* Model Version Parameters */
	INT4  IMRPhenomXInspiralPhaseVersion;
	INT4  IMRPhenomXIntermediatePhaseVersion;
	INT4  IMRPhenomXRingdownPhaseVersion;

	INT4  IMRPhenomXInspiralAmpVersion;
	INT4  IMRPhenomXIntermediateAmpVersion;
	INT4  IMRPhenomXRingdownAmpVersion;

	/* Mass Parameters */
	REAL8 m1_SI; 		// Mass in SI units
	REAL8 m2_SI;	 	// Mass in SI units
	REAL8 Mtot_SI; 	// Total mass in SI units
	REAL8 m1;		 		// Mass in solar masses
	REAL8 m2;	   		// Mass in solar masses
	REAL8 Mtot;  		// Total mass in solar masses
	REAL8 Mc;				// Chirp mass in solar masses
	REAL8 q;				// Mass ratio >= 1
	REAL8 eta;			// Symmetric mass ratio
	REAL8 delta;		// PN symmetry parameter: sqrt(1-4*eta)

	/* Spin Parameters */
	REAL8 chi1L;
	REAL8 chi2L;
	REAL8 chiEff;
	REAL8 chiPNHat;
	REAL8 STotR;
	REAL8 dchi;
	REAL8 SL;
	REAL8 SigmaL;

	/* Useful Powers (?) */
	REAL8 eta2;
	REAL8 eta3;
	REAL8 eta4;
	REAL8 chi1L2;
	REAL8 chi1L3;
	REAL8 chi2L2;
	REAL8 chi2L3;
	REAL8 chi1L2L; // chi1 * chi2;

	/* Amplitude and Phase Normalisation */
	REAL8 dphase0;
	REAL8 amp0;
	REAL8 ampNorm;

	/* Frequencies */
	REAL8 fMECO;
	REAL8 fISCO;

	/* Ringdown and Damping Frequencies for 22 Mode */
	REAL8 fRING;
	REAL8 fDAMP;

	/* Ringdown and Damping Frequencies for 21 Mode */
	REAL8 fRING21;
	REAL8 fDAMP21;

	/* Ringdown and Damping Frequencies for 33 Mode */
	REAL8 fRING33;
	REAL8 fDAMP33;

	/* Ringdown and Damping Frequencies for 32 Mode */
	REAL8 fRING32;
	REAL8 fDAMP32;

	/* Ringdown and Damping Frequencies for 44 Mode */
	REAL8 fRING44;
	REAL8 fDAMP44;

	REAL8 fMin;
	REAL8 fMax;
	REAL8 f_max_prime;
	REAL8 deltaF;
	REAL8 fCut;

	// Dimensionless frequency (Mf) defining the end of the waveform
	REAL8 fCutDef;

	REAL8 fRef;
	REAL8 MfRef;
	REAL8 M_sec;
	REAL8 phiRef_In;
	REAL8 phi0;
	REAL8 phifRef;
	REAL8 piM;
	REAL8 v_ref;

	/* Final mass and spin */
	REAL8 Erad;
	REAL8 afinal;
	REAL8 Mfinal;

	REAL8 distance;
	REAL8 inclination;
	REAL8 beta;

} IMRPhenomXWaveformStruct;


typedef struct tagIMRPhenomX_UsefulPowers
{
	REAL8 seven_sixths;
	REAL8 one_sixth;
	REAL8 eight_thirds;
	REAL8 seven_thirds;
	REAL8 five_thirds;
	REAL8 four_thirds;
	REAL8 two_thirds;
	REAL8 one_third;
	REAL8 five;
	REAL8 four;
	REAL8 three;
	REAL8 two;
	REAL8 sqrt;
	REAL8 itself;
	REAL8 m_sqrt;
	REAL8 m_one;
	REAL8 m_two;
	REAL8 m_three;
	REAL8 m_four;
	REAL8 m_five;
	REAL8 m_six;
	REAL8 m_one_third;
	REAL8 m_two_thirds;
	REAL8 m_four_thirds;
	REAL8 m_five_thirds;
	REAL8 m_seven_thirds;
	REAL8 m_eight_thirds;
	REAL8 m_one_sixth;
	REAL8 m_seven_sixths;

	REAL8 log;

} IMRPhenomX_UsefulPowers;

/*
 * useful powers of LAL_PI, calculated once and kept constant - to be initied with a call to
 */
extern IMRPhenomX_UsefulPowers powers_of_lalpi;

typedef struct tagIMRPhenomXPhaseCoefficients
{
	/* PHASE */
	/* Phase Transition Frequencies */
	REAL8 fPhaseInsMin;
	REAL8 fPhaseInsMax;
	REAL8 fPhaseIntMin;
	REAL8 fPhaseIntMax;
	REAL8 fPhaseRDMin;
	REAL8 fPhaseRDMax;

	REAL8 fPhaseMatchIN;
	REAL8 fPhaseMatchIM;

	REAL8 C1Int, C2Int;
	REAL8 C1MRD, C2MRD;

  /* These are the RD phenomenological coefficients 					*/
  REAL8 c0, c1, c2, c3, c4, cL, cRD;

  /* These are the intermediate phenomenological coefficients */
  REAL8 b0, b1, b2, b3, b4;

  /* These are the inspiral phenomenological coefficients 		*/
  REAL8 a0, a1, a2, a3, a4;

	/* Pre-cached variables */
	REAL8 c4ov3, cLovfda;

	/* TaylorF2 PN Coefficients */
	REAL8 phi0, phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10, phi11, phi12, phi13, phi5L, phi6L, phi8L, phi9L;
	REAL8 phi_initial, phiNorm;
	REAL8 dphi0, dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi8, dphi9, dphi10, dphi11, dphi12, dphi13, dphi5L, dphi6L, dphi8L, dphi9L;

	/* Pseudo-PN Coefficients */
	REAL8 sigma0, sigma1, sigma2, sigma3, sigma4, sigma5;

  /* Flag to set how many collocation points the RD region uses 	*/
  INT4  NCollocationPointsRD;

  /* Flag to set how many collocation points the INT region uses 	*/
  INT4  NCollocationPointsInt;

  /* Integer to tell us how many pseudo PN terms are used 											*/
  INT4	NPseudoPN;
	INT4  NCollocationPointsPhaseIns;

	/* The canonical ringdown phase is constructed from 5 collocation points 			*/
	REAL8 CollocationPointsPhaseRD[N_MAX_COLLOCATION_POINTS_PHASE_RD];
	REAL8 CollocationValuesPhaseRD[N_MAX_COLLOCATION_POINTS_PHASE_RD];
	REAL8 CoefficientsPhaseRD[N_MAX_COLLOCATION_POINTS_PHASE_RD];


	/* The canonical intermediate phase is constructed from 4/5 collocation points  */
	REAL8 CollocationPointsPhaseInt[N_MAX_COLLOCATION_POINTS_PHASE_INT];
	REAL8 CollocationValuesPhaseInt[N_MAX_COLLOCATION_POINTS_PHASE_INT];
	REAL8 CoefficientsPhaseInt[N_MAX_COLLOCATION_POINTS_PHASE_INT];

	/*
		 	For N pseudo-PN terms we need N+1 collocation points:
		 	We have set N_MAX_COLLOCATION_POINTS_INS = 5 to allow
		 	either 3 or 4 pseudo-PN coefficients to be used.
	*/
	REAL8 CollocationPointsPhaseIns[N_MAX_COLLOCATION_POINTS_PHASE_INS];
	REAL8 CollocationValuesPhaseIns[N_MAX_COLLOCATION_POINTS_PHASE_INS];
	REAL8 CoefficientsPhaseIns[N_MAX_COLLOCATION_POINTS_PHASE_INS];

} IMRPhenomXPhaseCoefficients;

typedef struct tagIMRPhenomXAmpCoefficients
{
	REAL8 fAmpInsMin;
	REAL8 fAmpInsMax;
	REAL8 fAmpIntMin;
	REAL8 fAmpIntMax;
	REAL8 fAmpRDMin;
	REAL8 fAmpRDMax;

	REAL8 fAmpMatchIN;
	REAL8 fAmpMatchIM;

  /* These are the RD phenomenological coefficients 					*/
  REAL8 c0, c1, c2, c3, c4, cL;

  /* These are the intermediate phenomenological coefficients */
  REAL8 b0, b1, b2, b3, b4, b5;

  /* These are the inspiral phenomenological coefficients 		*/
  REAL8 a0, a1, a2, a3, a4, a5;

	REAL8 v1RD, sigmaRD;

	REAL8 rho1, rho2, rho3;																	/* Inspiral pseudo-PN coefficients 							*/
	REAL8 delta0, delta1, delta2, delta3, delta4, delta5;		/* Intermediate phenomenological coefficients 	*/
	REAL8 gamma1, gamma2, gamma3;														/* Ringdown phenomenological coefficients 			*/
	REAL8 gammaR, gammaD13, gammaD2;

	/* PN Amplitude Prefactors */
	REAL8 pnInitial, pnOneThird, pnTwoThirds, pnThreeThirds, pnFourThirds, pnFiveThirds, pnSixThirds, pnSevenThirds, pnEightThirds, pnNineThirds;

  /* Flags to set the ringdown amplitude version			*/
	INT4  NCollocationPointsRD;
	INT4  IMRPhenomXRingdownAmpVersion;

	/* Flags to set the intermediate amplitude version	*/
	INT4  NCollocationPointsInt;
	INT4  IMRPhenomXIntermediateAmpVersion;

  /* Flags to set the inspiral amplitude version 			*/
	INT4	NPseudoPN;
	INT4  IMRPhenomXInspiralAmpVersion;

	/* The ringdown is constructed from 5 collocation points 						*/
	REAL8 CollocationPointsAmpRD[N_MAX_COLLOCATION_POINTS_AMP_RD];
	REAL8 CollocationValuesAmpRD[N_MAX_COLLOCATION_POINTS_AMP_RD];

	REAL8 CollocationPointsAmpInt[N_MAX_COLLOCATION_POINTS_AMP_INT];
	REAL8 CollocationValuesAmpInt[N_MAX_COLLOCATION_POINTS_AMP_INT];

	/* For 3 pseudo-PN parameters we need 4 collocation points         	*/
	/* For 4 pseudo-PN parameters we need 5 collocation points, etc... 	*/
	REAL8 CollocationPointsAmpIns[N_MAX_COLLOCATION_POINTS_AMP_INS];
	REAL8 CollocationValuesAmpIns[N_MAX_COLLOCATION_POINTS_AMP_INS];

} IMRPhenomXAmpCoefficients;



///////////////////////////// Useful Numerical Routines /////////////////////////////
int IMRPhenomX_Initialize_Powers(IMRPhenomX_UsefulPowers *p, REAL8 number);
int IMRPhenomX_Initialize_Powers_Light(IMRPhenomX_UsefulPowers *p, REAL8 number);

int IMRPhenomXSetWaveformVariables(
IMRPhenomXWaveformStruct *pWF,
	const REAL8 m1_SI,
	const REAL8 m2_SI,
	const REAL8 chi1L_In,
	const REAL8 chi2L_In,
	const REAL8 deltaF,
	const REAL8 fRef,
	const REAL8 phiRef,
	const REAL8 f_min,
	const REAL8 f_max,
	const REAL8 distance,
	const REAL8 inclination,
	LALDict *lalParams,
	const UINT4 debug
);

double IMRPhenomX_Amplitude_22(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXAmpCoefficients *pAmp, IMRPhenomXWaveformStruct *pWF);
double IMRPhenomX_Phase_22(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *Phase, IMRPhenomXWaveformStruct *pWF);
double IMRPhenomX_dPhase_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF);

int IMRPhenomXGetAmplitudeCoefficients(IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp);
int IMRPhenomXGetPhaseCoefficients(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);

double IMRPhenomX_TimeShift_22(IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF);
    
void IMRPhenomX_Phase_22_ConnectionCoefficients(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase);

/* Function to check if the input mode array contains supported modes */
INT4 check_input_mode_array(LALDict *lalParams);

#ifdef __cplusplus
}
#endif

#endif	// of #ifndef _LALSIM_IMR_PHENOMX_INTERNALS_H
