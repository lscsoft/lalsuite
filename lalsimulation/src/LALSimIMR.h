/*
 * Copyright (C) 2011 N. Fotopoulos <nickolas.fotopoulos@ligo.org>
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

#ifndef _LALSIMIMR_H
#define _LALSIMIMR_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @defgroup LALSimIMR_h Header LALSimIMR.h
 * @ingroup lalsimulation_inspiral
 *
 * @brief Routines for generating inspiral-merger-ringdown waveforms.
 *
 * @{
 * @defgroup LALSimIMRPhenom_c                   LALSimIMRPhenom.c
 * @defgroup LALSimIMREOBNRv2_c                  LALSimIMREOBNRv2.c
 * @defgroup LALSimIMRSpinAlignedEOB_c           LALSimIMRSpinAlignedEOB.c
 * @defgroup LALSimIMRSEOBNRROM_c                LALSimIMRSEOBNRvxROMXXX.c
 * @defgroup LALSimIMRSEOBNRv2ChirpTime_c        LALSimIMRSEOBNRv2ChirpTime.c
 * @defgroup LALSimIMRPSpinInspiralRD_c          LALSimIMRPSpinInspiralRD.c
 * @}
 *
 * @addtogroup LALSimIMR_h
 * @{
 */

/**
 * The number of e-folds of ringdown which should be attached for
 * EOBNR models
 */
#define EOB_RD_EFOLDS 10.0

/** @} */

/* in module LALSimIMRPhenom.c */

int XLALSimIMRPhenomAGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phiPeak, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_min, const REAL8 f_max, const REAL8 distance);
int XLALSimIMRPhenomAGenerateTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const REAL8 inclination);
double XLALSimIMRPhenomBComputeChi(const REAL8 m1, const REAL8 m2, const REAL8 s1z, const REAL8 s2z);
double XLALSimIMRPhenomAGetFinalFreq(const REAL8 m1, const REAL8 m2);
double XLALSimIMRPhenomBGetFinalFreq(const REAL8 m1, const REAL8 m2, const REAL8 chi);
int XLALSimIMRPhenomBGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phiPeak, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance);
int XLALSimIMRPhenomBGenerateTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const REAL8 inclination);
int XLALSimIMRPhenomBMetricInMEtaChi(REAL8 *gamma00, REAL8 *gamma01, REAL8 *gamma02, REAL8 *gamma11, REAL8 *gamma12, REAL8 *gamma22, const REAL8 m1, const REAL8 m2, const REAL8 chi, const REAL8 fLow, const REAL8FrequencySeries *Sh);
int XLALSimIMRPhenomBMetricInTheta0Theta3Theta3S(REAL8 *gamma00, REAL8 *gamma01, REAL8 *gamma02, REAL8 *gamma11, REAL8 *gamma12, REAL8 *gamma22, const REAL8 m1, const REAL8 m2, const REAL8 chi, const REAL8 fLow, const REAL8FrequencySeries *Sh);
double XLALSimIMRPhenomCGetFinalFreq(const REAL8 m1, const REAL8 m2, const REAL8 chi);
int XLALSimIMRPhenomCGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phiPeak, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance);
int XLALSimIMRPhenomCGenerateTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const REAL8 inclination);
int XLALSimIMRPhenomDGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi1, const REAL8 chi2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance);
int XLALSimIMRPhenomP(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 chi_eff, const REAL8 chip, const REAL8 eta, const REAL8 thetaJ, const REAL8 Mtot_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, const REAL8 f_ref);
int XLALSimIMRPhenomPFrequencySequence(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, const REAL8 chi_eff, const REAL8 chip, const REAL8 eta, const REAL8 thetaJ, const REAL8 Mtot_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 f_ref);
int XLALSimIMRPhenomPCalculateModelParameters(REAL8 *chi_eff, REAL8 *chip, REAL8 *eta, REAL8 *thetaJ, REAL8 *alpha0, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 lnhatx, const REAL8 lnhaty, const REAL8 lnhatz, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z);

/* in module LALSimIMREOBNRv2.c */

int XLALSimIMREOBNRv2DominantMode(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
int XLALSimIMREOBNRv2AllModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
SphHarmTimeSeries *XLALSimIMREOBNRv2Modes(const REAL8 phiRef, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 fLower, const REAL8 distance);


/* in module LALSimIMRSpinAlignedEOB.c */

double XLALSimIMRSpinAlignedEOBPeakFrequency(REAL8 m1SI, REAL8 m2SI, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion);
int XLALSimIMRSpinAlignedEOBWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 inc, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion);
//int XLALSimIMRSpinEOBWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 inc, const REAL8 spin1[], const REAL8 spin2[]);


/* in module LALSimIMRSEOBNRv1ROMEffectiveSpin.c */

int XLALSimIMRSEOBNRv1ROMEffectiveSpin(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi);
int XLALSimIMRSEOBNRv1ROMEffectiveSpinFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi);

/* in module LALSimIMRSEOBNRv1ROMDoubleSpin.c */

int XLALSimIMRSEOBNRv1ROMDoubleSpin(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv1ROMDoubleSpinFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);


/* in module LALSimIMRSEOBNRv2ROMEffectiveSpin.c */

int XLALSimIMRSEOBNRv2ROMEffectiveSpin(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi);
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi);
int XLALSimIMRSEOBNRv2ROMEffectiveSpinTimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi);
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi);


/* in module LALSimIMRSEOBNRv2ROMDoubleSpin.c */

int XLALSimIMRSEOBNRv2ROMDoubleSpin(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv2ROMDoubleSpinTimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);


/* in module LALSimIMRSEOBNRv2ChirpTime.c */

REAL8 XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min);


/* in module LALSimIMRPSpinInspiralRD.c */

int XLALSimIMRPhenSpinFinalMassSpin(REAL8 *finalMass, REAL8 *finalSpin, REAL8 m1, REAL8 m2, REAL8 s1s1, REAL8 s2s2, REAL8 s1L, REAL8 s2L, REAL8 s1s2, REAL8 energy); 
int XLALSimSpinInspiralGenerator(REAL8TimeSeries **hPlus, REAL8TimeSeries **hCross, REAL8 phi_start, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *testGRparams);
int XLALSimIMRPhenSpinInspiralRDGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phi0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralWaveformFlags *waveFlag, LALSimInspiralTestGRParam *testGRparam);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMR_H */
