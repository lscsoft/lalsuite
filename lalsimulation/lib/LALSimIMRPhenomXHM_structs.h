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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

//
//  LALSimIMRPhenomXHM_structs.h
//
//  Created by Marta on 13/02/2019.
//
#ifndef LALSimIMRPhenomXHM_structs_h
#define LALSimIMRPhenomXHM_structs_h

#ifdef __cplusplus
extern "C" {
#endif


#include <lal/LALAtomicDatatypes.h>
#include "LALSimIMRPhenomX_internals.h"

#define N_HIGHERMODES_IMPLEMENTED 4 // lm = (21, 33, 32, 44)

#define N_MAX_COEFFICIENTS_PHASE_INS 15       //Maximun number of coefficients of the inspiral ansatz
#define N_MAX_COEFFICIENTS_PHASE_INTER 6      //Maximun number of coefficients of the intermediate ansatz
#define N_MAX_COEFFICIENTS_PHASE_RING 8       //Maximun number of coefficients of the ringdown ansatz. 5 + 3 optional for flattened phi'

#define N_MAX_COEFFICIENTS_AMPLITUDE_INS 3    //Maximun number of collocation points in the inspiral
#define N_MAX_COEFFICIENTS_AMPLITUDE_INTER 8  //Maximun number of collocation points in the intermediate. New release 4 coll points + 2x2 boundaries
#define N_MAX_COEFFICIENTS_AMPLITUDE_RING 6   //Maximun number of coefficients in the ringdown. We only use 3 degrees of freedom, but we have the double to store fits of 3 coefficients or 3 collocation points
#define N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX 4  //Maximun number of coefficients in the ringdown auxiliar region for mode-mixing. Default 2 collocation points + point & derivative right boundary

// Data structure to hold QNM frequencies
typedef double (*fitQNM_fring) (double finalDimlessSpin);
typedef double (*fitQNM_fdamp) (double finalDimlessSpin);

typedef struct tagQNMFits {
 	fitQNM_fring fring_lm[N_HIGHERMODES_IMPLEMENTED];
 	fitQNM_fdamp fdamp_lm[N_HIGHERMODES_IMPLEMENTED];
} QNMFits;


// General fit function. This is a type for defining the functions for the parameter space fits.
typedef double (*ParameterSpaceFit) (IMRPhenomXWaveformStruct *pWF, int flag);

// Waveform struct.  Store useful variable specific of the higher modes that are not in 22 IMRPhenomXWaveformStruct or that need to be updated mode-by-mode. */
typedef struct tagIMRPhenomXHMWaveformStruct
{
        /* Model Version Parameters */
        INT4  IMRPhenomXHMInspiralAmpVersion;
        INT4  IMRPhenomXHMIntermediateAmpVersion;
        INT4  IMRPhenomXHMRingdownAmpVersion;

        INT4  IMRPhenomXHMInspiralPhaseVersion;
        INT4  IMRPhenomXHMIntermediatePhaseVersion;
        INT4  IMRPhenomXHMRingdownPhaseVersion;

        INT4  IMRPhenomXHMInspiralAmpFitsVersion;
        INT4  IMRPhenomXHMIntermediateAmpFitsVersion;
        INT4  IMRPhenomXHMRingdownAmpFitsVersion;

        INT4  IMRPhenomXHMInspiralPhaseFitsVersion;
        INT4  IMRPhenomXHMIntermediatePhaseFitsVersion;
        INT4  IMRPhenomXHMRingdownPhaseFitsVersion;

        INT4  IMRPhenomXHMInspiralAmpFreqsVersion;
        INT4  IMRPhenomXHMIntermediateAmpFreqsVersion;
        INT4  IMRPhenomXHMRingdownAmpFreqsVersion;

        INT4  IMRPhenomXHMInspiralPhaseFreqsVersion;
        INT4  IMRPhenomXHMIntermediatePhaseFreqsVersion;
        INT4  IMRPhenomXHMRingdownPhaseFreqsVersion;



        INT4 IMRPhenomXHMReleaseVersion;

        /* Parameters that define deviation of the tuned coprecessing mode PhenomXCP from PhenomX (500) */
        REAL8 MU1;   // MR Amplitude
        REAL8 MU2;   // MR Amplitude: modifies gamma2
        REAL8 MU3;   // MR Amplitude: modifies gamma3
        REAL8 MU4;   // MR Amplitude: would modify appearance of fRing in MR amplitude
        REAL8 NU0;   // MR Phase
        REAL8 NU4;   // MR Phase
        REAL8 NU5;   // MR Phase
        REAL8 NU6;   // MR Phase
        REAL8 ZETA1; // INT Phase
        REAL8 ZETA2; // INT Phase
        REAL8 PNR_DEV_PARAMETER; // Is zero when no precession, and non-zero otherwise


        /* Spin Parameters */
        REAL8 chi_s, chi_a;  // (chi1 +/- chi2)/2

        /* MECO, Ringdown and Damping Frequencies */
		    REAL8 fMECOlm;  // = wf22->fMECO*m/2
        REAL8 fRING;
        REAL8 fDAMP;

        /* Limit between comparable and extreme mass ratios for the phase */
        REAL8 etaEMR;

        /* The PhenomXAS functions returns the phase without a linear part that is not relevant for the matches.
          We need this linear part for the multimode waveform. timeshift * f + phaseshift */
        REAL8 timeshift;
        REAL8 phaseshift;
        /* NOTE that we wish to destinguish between the global timshihft for 22, and that for each higher moment */
        REAL8 timeshiftLM;
        REAL8 phaseshiftLM;
        REAL8 phiref22;  //Correction to apply to PhX phase for returning phiRef at reference frequency

        /* mode labels ad tags*/
        INT4 ell;      // spherical harmomic l
        INT4 emm;      // spherical harmomic m
        INT4 modeInt;
        INT4 modeTag;

        /* Number collocation points */
        INT4 nCollocPtsInterPhase;
        INT4 nCollocPtsRDPhase;
        INT4 nCollocPtsInspAmp;
        INT4 nCollocPtsInterAmp;

        /* variable to control use of mode-mixing routines */
        INT4 MixingOn;

        // 2x2 matrix containing the mixing coefficients of the lm mode with the nearest dominant neighbour
        COMPLEX8 mixingCoeffs[4];
        COMPLEX8 mu22, mu32;

        /* Variables to control use of Amplitude Veto for removing certain collocation points. */
        INT4 InspiralAmpVeto, IntermediateAmpVeto, RingdownAmpVeto;

        // The same than in WF22, sqrt(2*eta/3)/Pi, this is the normalization factor of the leading order of the 22.
        REAL8 ampNorm;

        // Variable to control if the case is EMR and apply the two intermediate regions for the amplitude
        INT4 AmpEMR;

        // Variable to control if the odd modes are zero
        INT4 Ampzero;

        // Multiply by Amp0 to transform NR amplitude to physical amplitude (= pWF22->ampNorm * pWF22->amp0 )
        REAL8 Amp0;

        // Variable to control the use of FAmpPN function for 21 instead of the power series
        INT4 useFAmpPN;

        /* time-shift of the peak of the hybrids' 22 wrt end of the waveform*/
        REAL8 DeltaT;

        REAL8 fPhaseRDflat;      // Ringdown -> Flattened region
        REAL8 fAmpRDfalloff;     // Ringdown -> Falloff region

    } IMRPhenomXHMWaveformStruct;


    /******************* AMPLITUDE ****************/

    // Store the coefficients needed to perform the amplitude reconstruction
    typedef struct tagIMRPhenomXHMAmpCoefficients
    {
            /* Cutting frequencies */
            REAL8 fAmpMatchIN;    // Inspiral -> Intermediate
            REAL8 fAmpMatchInt12; // Intermediate1 -> Intermediate2. Only for EMR cases
            REAL8 fAmpMatchIM;    // Intermediate -> Ringdown
            REAL8 fRDAux;  // Auxiliar Ringdown region for mode-mixing

            /* PN Amplitude Prefactors */
            COMPLEX16 pnInitial, pnOneThird, pnTwoThirds, pnThreeThirds, pnFourThirds, pnFiveThirds, pnSixThirds, pnSevenThirds, pnEightThirds,pnNineThirds;

            /* PN Amplitude dominant factor = Pi * Sqrt(2 eta/3) (2Pi /m)^(-7/6) */
            REAL8 PNdominant, PNdominantlm, ampNorm;

            UINT2 PNdominantlmpower;

            /* PN Amplitude global prefactor = leading lm order after substracting PNdominant*/
            REAL8 PNglobalfactor;

            /* Coefficients of the pseudo-PN terms */
            REAL8 rho1, rho2, rho3;

            /* Coefficients of the polynomial in the intermediate region. 5th order by default */
            REAL8 delta0, delta1, delta2, delta3, delta4, delta5;

            /* Coefficients of the polynomial in the first intermediate region (for EMR only). 4th order.*/
            REAL8 alpha0, alpha1, alpha2, alpha3, alpha4;

            // fits of coefficients/collocation points
            ParameterSpaceFit InspiralAmpFits[N_HIGHERMODES_IMPLEMENTED*N_MAX_COEFFICIENTS_AMPLITUDE_INS];
            ParameterSpaceFit IntermediateAmpFits[N_HIGHERMODES_IMPLEMENTED*N_MAX_COEFFICIENTS_AMPLITUDE_INTER];
            ParameterSpaceFit RingdownAmpFits[N_HIGHERMODES_IMPLEMENTED*N_MAX_COEFFICIENTS_AMPLITUDE_RING + N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX];

            /* Flag to set how many collocation points the inspiral region uses  */
            REAL8 CollocationPointsValuesAmplitudeInsp[N_MAX_COEFFICIENTS_AMPLITUDE_INS];
            REAL8 CollocationPointsFreqsAmplitudeInsp[N_MAX_COEFFICIENTS_AMPLITUDE_INS];
            REAL8 InspiralCoefficient[N_MAX_COEFFICIENTS_AMPLITUDE_INS];

            /* Flag to set how many collocation points the intermediate region uses */
            REAL8 CollocationPointsFreqsAmplitudeInter[N_MAX_COEFFICIENTS_AMPLITUDE_INTER];
            REAL8 CollocationPointsValuesAmplitudeInter[N_MAX_COEFFICIENTS_AMPLITUDE_INTER];
            REAL8 InterCoefficient[N_MAX_COEFFICIENTS_AMPLITUDE_INTER];
            UINT2 nCoefficientsInter;
            UINT2 VersionCollocPtsInter[N_MAX_COEFFICIENTS_AMPLITUDE_INTER];

            /* Flag to set how many collocation points the intermediate region uses */
            REAL8 CollocationPointsFreqsAmplitudeRD[N_MAX_COEFFICIENTS_AMPLITUDE_RING];
            REAL8 CollocationPointsValuesAmplitudeRD[N_MAX_COEFFICIENTS_AMPLITUDE_RING];
            REAL8 RDCoefficient[N_MAX_COEFFICIENTS_AMPLITUDE_RING];
            REAL8 CollocationPointsFreqsAmplitudeRDAux[N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX];
            REAL8 CollocationPointsValuesAmplitudeRDAux[N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX];
            REAL8 RDAuxCoefficient[N_MAX_COEFFICIENTS_AMPLITUDE_RDAUX];
            UINT2 nCollocPtsRDAux, nCoefficientsRDAux;

            // Frequencies, values and derivatives for the intermediate reconstruction
            // The frequencies are the same than in CollocationPointsFreqsAmplitudeInsp for the corresponding mode,
            // but in this way they are easier accesible.
            REAL8 f1,f2,f3,f4,v1,v2,v3,v4,d1,d4;

            // Order of the polynomial in the intermediate region. 5th->105, for the first EMR region is 1042
            INT4 InterAmpPolOrder;

        		// Store the PN amplitude at the frequencies of the collocation points in the inspiral
        		REAL8 PNAmplitudeInsp[N_MAX_COEFFICIENTS_AMPLITUDE_INS];

            // For the pseudo part of Inspiral Amplitude ansatz. Used in LALSimIMRPhenomXHM_inspiral.c
            REAL8 fcutInsp_seven_thirds;
            REAL8 fcutInsp_eight_thirds;
            REAL8 fcutInsp_three;

            // Coefficients of the phasing factor coming from the SPA and time-domain Post-newtonian series, only for the 21 mode
            REAL8 xdot5, xdot6, xdot65, xdot7, xdot75, xdot8, xdot8Log, xdot85, log2pi_two_thirds;
            COMPLEX16 x05, x1, x15, x25, x2, x3;
            REAL8 PNTDfactor;

            // Variable to control the use of the inspiral ansatz to compute the amplitude at fmaxCalc
            INT4 useInspAnsatzRingdown;

            // Variables to control if we have to check that the collocation points are wavy
            INT4 WavyInsp, WavyInt;

            // Amp0 = wf22->ampNorm * wf22->amp0. Multiplying by this gives the amp factor of the 22 and transforms to "physical" units
            REAL8 Amp0;

            UINT2 InspRescaleFactor, InterRescaleFactor, RDRescaleFactor;

    } IMRPhenomXHMAmpCoefficients;


/*************** PHASE *****************/

// Store the coefficients needed to perform the phase reconstruction
typedef struct tagIMRPhenomXHMPhaseCoefficients
{
        /* Phase Transition Frequencies */
        REAL8 fPhaseMatchIN;
        REAL8 fPhaseMatchIM;

        REAL8 deltaphiLM;

        /* These are the RD phenomenological coefficients, with mode-mixing off */
        REAL8 alpha0, alpha2, alphaL;
        REAL8 phi0RD, dphi0RD;

        /* These are the RD phenomenological coefficients, with mode-mixing on: they give a phenomenological representation of the spheroidal-harmonic phase */
        REAL8 phi0_S, alpha0_S, alpha2_S, alpha4_S, alphaL_S;

        /* These are the intermediate phenomenological coefficients */
        REAL8 c0, c1, c2, c3, c4, cL;
        REAL8 CRD, C1RD, CINSP, C1INSP;

        /* These are the inspiral phenomenological coefficients */
        REAL8 phi[N_MAX_COEFFICIENTS_PHASE_INS];
        REAL8 phiL[N_MAX_COEFFICIENTS_PHASE_INS];
        REAL8 LambdaPN;

        // fits of coefficients/collocation points
        ParameterSpaceFit InspiralPhaseFits[N_HIGHERMODES_IMPLEMENTED];
        ParameterSpaceFit IntermediatePhaseFits[N_HIGHERMODES_IMPLEMENTED*N_MAX_COEFFICIENTS_PHASE_INTER];
        ParameterSpaceFit RingdownPhaseFits[N_MAX_COEFFICIENTS_PHASE_RING];

        /* Flag to set how many collocation points the RD region uses --used only for the 32 mode */
        INT4  NCollocationPointsRD;
        REAL8 CollocationPointsValuesPhaseRD[N_MAX_COEFFICIENTS_PHASE_RING];
        REAL8 CollocationPointsFreqsPhaseRD[N_MAX_COEFFICIENTS_PHASE_RING];
        REAL8 RDCoefficient[N_MAX_COEFFICIENTS_PHASE_RING + 3]; // Add the three coefficient for rational decay

        /* Flag to set how many collocation points the intermediate region uses */
        INT4  NCollocationPointsInt;
        REAL8 CollocationPointsFreqsPhaseInter[N_MAX_COEFFICIENTS_PHASE_INTER];
        REAL8 CollocationPointsValuesPhaseInter[N_MAX_COEFFICIENTS_PHASE_INTER];

} IMRPhenomXHMPhaseCoefficients;


#ifdef __cplusplus
}
#endif

#endif /* LALSimIMRPhenomXHM_structs_h */
