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

#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
#endif

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
 * @defgroup LALSimIMRSpinPrecEOB_c              LALSimIMRSpinPrecEOB.c
 * @defgroup LALSimIMRSpinPrecEOBv4P_c           LALSimIMRSpinPrecEOBv4P.c
 * @defgroup LALSimIMRSEOBNRROM_c                LALSimIMRSEOBNRvxROMXXX.c
 * @defgroup LALSimIMRSEOBNRv2ChirpTime_c        LALSimIMRSEOBNRv2ChirpTime.c
 * @defgroup LALSimIMRPSpinInspiralRD_c          LALSimIMRPSpinInspiralRD.c
 * @defgroup LALSimIMRTidal_c                    LALSimIMRLackeyTidal2013.c
 * @defgroup LALSimNRSur7dq2_c                   LALSimIMRNRSur7dq2.c
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

typedef enum tagIMRPhenomP_version_type {
 IMRPhenomPv1_V, /**< version 1: based on IMRPhenomC */
 IMRPhenomPv2_V,  /**< version 2: based on IMRPhenomD */
 IMRPhenomPv2NRTidal_V /**< version Pv2_NRTidal: based on IMRPhenomPv2; NRTides (https://arxiv.org/pdf/1706.02969.pdf) added before precession */
} IMRPhenomP_version_type;

typedef enum tagSEOBNRv4TSurrogate_spline_order {
  SEOBNRv4TSurrogate_CUBIC, /**< use cubic splines in frequency */
  SEOBNRv4TSurrogate_LINEAR /**< use linear splines in frequency */
} SEOBNRv4TSurrogate_spline_order;

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
int XLALSimIMRPhenomCGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phiPeak, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, LALDict *extraParams);
int XLALSimIMRPhenomCGenerateTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const REAL8 inclination, LALDict *extraParams);

/* in module LALSimIMRPhenomD.c */
int XLALSimIMRPhenomDGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 fRef, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi1, const REAL8 chi2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, LALDict *extraParams);
int XLALSimIMRPhenomDFrequencySequence(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi0, const REAL8 fRef_in, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi1, const REAL8 chi2, const REAL8 distance, LALDict *extraParams);
double XLALIMRPhenomDGetPeakFreq(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in);
double XLALSimIMRPhenomDChirpTime(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in, const REAL8 fHz);
double XLALSimIMRPhenomDFinalSpin(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in);

int XLALSimIMRPhenomP(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 chi1_l, const REAL8 chi2_l, const REAL8 chip, const REAL8 thetaJ, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, const REAL8 f_ref, IMRPhenomP_version_type IMRPhenomP_version, LALDict *extraParams);
int XLALSimIMRPhenomPFrequencySequence(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, const REAL8 chi1_l, const REAL8 chi2_l, const REAL8 chip, const REAL8 thetaJ, REAL8 m1_SI, const REAL8 m2_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 f_ref, IMRPhenomP_version_type IMRPhenomP_version, LALDict *extraParams);
int XLALSimIMRPhenomPCalculateModelParametersOld(REAL8 *chi1_l, REAL8 *chi2_l, REAL8 *chip, REAL8 *thetaJ, REAL8 *alpha0, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 lnhatx, const REAL8 lnhaty, const REAL8 lnhatz, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, IMRPhenomP_version_type IMRPhenomP_version);
int XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(REAL8 *chi1_l, REAL8 *chi2_l, REAL8 *chip, REAL8 *thetaJN, REAL8 *alpha0, REAL8 *phi_aligned, REAL8 *zeta_polariz, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 phiRef, const REAL8 incl, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, IMRPhenomP_version_type IMRPhenomP_version);

/* in module LALSimIMREOBNRv2.c */

int XLALSimIMREOBNRv2DominantMode(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
int XLALSimIMREOBNRv2AllModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
SphHarmTimeSeries *XLALSimIMREOBNRv2Modes(const REAL8 phiRef, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 fLower, const REAL8 distance);


/* in module LALSimIMRSpinAlignedEOB.c */

double XLALSimIMRSpinAlignedEOBPeakFrequency(REAL8 m1SI, REAL8 m2SI, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion);
int XLALSimIMRSpinAlignedEOBWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 inc, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion, LALDict *LALparams);
int XLALSimIMRSpinAlignedEOBWaveformAll(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 inc, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion, const REAL8 lambda2Tidal1, const REAL8 lambda2Tidal2, const REAL8 omega02Tidal1, const REAL8 omega02Tidal2, const REAL8 lambda3Tidal1, const REAL8 lambda3Tidal2, const REAL8 omega03Tidal1, const REAL8 omega03Tidal2, const REAL8 quadparam1, const REAL8 quadparam2, REAL8Vector *nqcCoeffsInput, const INT4 nqcFlag, LALValue *ModeArray);
int XLALSimIMRSpinAlignedEOBModes(SphHarmTimeSeries ** hlmmode,
  //SM
  REAL8Vector ** dynamics_out, /**<< OUTPUT, low-sampling dynamics */
  REAL8Vector ** dynamicsHi_out, /**<< OUTPUT, high-sampling dynamics */
  //SM
 REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 spin1z, const REAL8 spin2z, UINT4 SpinAlignedEOBversion, const REAL8 lambda2Tidal1, const REAL8 lambda2Tidal2, const REAL8 omega02Tidal1, const REAL8 omega02Tidal2, const REAL8 lambda3Tidal1, const REAL8 lambda3Tidal2, const REAL8 omega03Tidal1, const REAL8 omega03Tidal2, const REAL8 quadparam1, const REAL8 quadparam2, REAL8Vector *nqcCoeffsInput, const INT4 nqcFlag);
/*int XLALSimIMRSpinEOBWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fMin, const REAL8 r, const REAL8 inc, const REAL8 spin1[], const REAL8 spin2[]);
 */

/* in module LALSimIMRSpinPrecEOB.c */

int XLALSimIMRSpinEOBWaveform(
                              REAL8TimeSeries **hplus,
                              REAL8TimeSeries **hcross,
                              const REAL8     phiC,
                              const REAL8     deltaT,
                              const REAL8     m1SI,
                              const REAL8     m2SI,
                              const REAL8     fMin,
                              const REAL8     r,
                              const REAL8     inc,
                              const REAL8     spin1[],
                              const REAL8     spin2[],
                              const UINT4     PrecEOBversion
                              );
int XLALSimIMRSpinEOBWaveformAll(
                                 REAL8TimeSeries **hplus,
                                 REAL8TimeSeries **hcross,
                                 REAL8Vector     **dynamicsHi,
                                 SphHarmTimeSeries **hlmPTSout,
                                 SphHarmTimeSeries **hlmPTSHi,
                                 SphHarmTimeSeries **hIMRlmJTSHi,
                                 SphHarmTimeSeries **hIMRoutput,
                                 REAL8Vector     **AttachParams,
                                 const REAL8     phiC,
                                 const REAL8     deltaT,
                                 const REAL8     m1SI,
                                 const REAL8     m2SI,
                                 const REAL8     fMin,
                                 const REAL8     r,
                                 const REAL8     inc,
                                 const REAL8     INspin1x,
                                 const REAL8     INspin1y,
                                 const REAL8     INspin1z,
                                 const REAL8     INspin2x,
                                 const REAL8     INspin2y,
                                 const REAL8     INspin2z,
                                 const UINT4     PrecEOBversion
                                 );
typedef enum tagflagSEOBNRv4P_hamiltonian_derivative {
  FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_ANALYTICAL = 0,  /**< use analytical derivatives (opt) */
  FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL  = 1  /**< use numerical derivatives (pre-opt) */
} flagSEOBNRv4P_hamiltonian_derivative;

typedef enum tagflagSEOBNRv4P_euler_extension {
  FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION = 0, /**< QNM-based simple precession prescription post-merger */
  FLAG_SEOBNRv4P_EULEREXT_CONSTANT = 1 /**< Euler angles set to constants post-merger */
} flagSEOBNRv4P_euler_extension;

typedef enum tagflagSEOBNRv4P_Zframe {
  FLAG_SEOBNRv4P_ZFRAME_L = 0, /**< set Z axis of the P-frame along L */
  FLAG_SEOBNRv4P_ZFRAME_LN = 1 /**< set Z axis of the P-frame along LN */
} flagSEOBNRv4P_Zframe;


/* in module LALSimIMRSpinPrecEOBv4P.c */

int XLALEOBHighestInitialFreq(REAL8 *freqMinRad, REAL8 mTotal);
int XLALEOBCheckNyquistFrequency(REAL8 m1, REAL8 m2, REAL8 spin1[3], REAL8 spin2[3],
			      UINT4 ell_max, Approximant approx,
			      REAL8 deltaT);
int XLALSimIMRSpinPrecEOBWaveformAll(
                                 REAL8TimeSeries   **hplus,
                                 REAL8TimeSeries   **hcross,
                                 SphHarmTimeSeries **hIlm,
                                 SphHarmTimeSeries **hJlm,
                                 REAL8Vector       **seobdynamicsAdaSVector,
                                 REAL8Vector       **seobdynamicsHiSVector,
                                 REAL8Vector       **seobdynamicsAdaSHiSVector,
                                 REAL8Vector       **tVecPmodes,
                                 REAL8Vector       **hP22_amp,
                                 REAL8Vector       **hP22_phase,
                                 REAL8Vector       **hP21_amp,
                                 REAL8Vector       **hP21_phase,
                                 REAL8Vector       **hP33_amp,
                                 REAL8Vector       **hP33_phase,
                                 REAL8Vector       **hP44_amp,
                                 REAL8Vector       **hP44_phase,
                                 REAL8Vector       **hP55_amp,
                                 REAL8Vector       **hP55_phase,
                                 REAL8Vector       **alphaJ2P,
                                 REAL8Vector       **betaJ2P,
                                 REAL8Vector       **gammaJ2P,
                                 REAL8Vector       **mergerParams,
                                 const REAL8       phiC,
                                 const REAL8       INdeltaT,
                                 const REAL8       m1SI,
                                 const REAL8       m2SI,
                                 const REAL8       fMin,
                                 const REAL8       r,
                                 const REAL8       inc,
                                 const REAL8       chi1x,
                                 const REAL8       chi1y,
                                 const REAL8       chi1z,
                                 const REAL8       chi2x,
                                 const REAL8       chi2y,
                                 const REAL8       chi2z,
                                 LALValue         *modearray,
                                 LALDict          *seobflags
                                 );

int XLALSimIMRSpinPrecEOBWaveform(
                              REAL8TimeSeries **hplus,
                              REAL8TimeSeries **hcross,
                              const REAL8     phiC,
                              const REAL8     deltaT,
                              const REAL8     m1SI,
                              const REAL8     m2SI,
                              const REAL8     fMin,
                              const REAL8     r,
                              const REAL8     inc,
                              const REAL8     INspin1[],
                              const REAL8     INspin2[],
                              const UINT4     PrecEOBversion,
                              LALDict         *LALParams
                              );

SphHarmTimeSeries *XLALSimIMRSpinPrecEOBModes(
                              const REAL8     deltaT,
                              const REAL8     m1SI,
                              const REAL8     m2SI,
                              const REAL8     fMin,
                              const REAL8     r,
                              const REAL8     INspin1[],
                              const REAL8     INspin2[],
                              const UINT4     PrecEOBversion,
                              LALDict         *LALParams
                              );

/* in module LALSimIMREOBNRv2HMROM.c */

int XLALSimIMREOBNRv2HMROM(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI,  REAL8 m2SI, const int higherModesFlag);

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
int XLALSimIMRSEOBNRv2ROMDoubleSpinAmpPhaseInterpolants( struct tagREAL8Vector **amplitude_interp, struct tagREAL8Vector **amplitude_freq_points, struct tagREAL8Vector **phase_interp, struct tagREAL8Vector **phase_freq_points, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);


/* in module LALSimIMRSEOBNRv2ROMDoubleSpinHI.c */

int XLALSimIMRSEOBNRv2ROMDoubleSpinHI(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max);
int XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max);
int XLALSimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);


/* in module LALSimIMRSEOBNRv2ChirpTime.c */

REAL8 XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 f_min);


/* in module LALSimIMRLackeyTidal2013.c */

int XLALSimIMRLackeyTidal2013(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 mBH_SI, REAL8 mNS_SI, REAL8 chi_BH, REAL8 Lambda);
int XLALSimIMRLackeyTidal2013FrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 mBH_SI, REAL8 mNS_SI, REAL8 chi_BH, REAL8 Lambda);


/* in module LALSimIMRSEOBNRv4ROM.c */

int XLALSimIMRSEOBNRv4ROM(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max);
int XLALSimIMRSEOBNRv4ROMFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max);
int XLALSimIMRSEOBNRv4ROMTimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv4ROMFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);

/* in module LALSimIMRSEOBNRv4ROM_NRTidal.c */

int XLALSimIMRSEOBNRv4ROMNRTidalFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 Lambda1, REAL8 Lambda2);
int XLALSimIMRSEOBNRv4ROMNRTidal(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 Lambda1, REAL8 Lambda2);

/* in module LALSimIMRSEOBNRv4TSurrogate.c */

int XLALSimIMRSEOBNRv4TSurrogate(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, SEOBNRv4TSurrogate_spline_order spline_order);
int XLALSimIMRSEOBNRv4TSurrogateFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, SEOBNRv4TSurrogate_spline_order spline_order);


/* in module LALSimIMRPSpinInspiralRD.c */

int XLALSimIMRPhenSpinFinalMassSpin(REAL8 *finalMass, REAL8 *finalSpin, REAL8 m1, REAL8 m2, REAL8 s1s1, REAL8 s2s2, REAL8 s1L, REAL8 s2L, REAL8 s1s2, REAL8 energy);
int XLALSimSpinInspiralGenerator(REAL8TimeSeries **hPlus, REAL8TimeSeries **hCross, REAL8 phi_start, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams);
int XLALSimIMRPhenSpinInspiralRDGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phi0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams);

/* in module LALSimInspiralNRWaveforms.c */

int XLALSimInspiralNRWaveformGetSpinsFromHDF5File(
  REAL8 *S1x,             /**< [out] Dimensionless spin1x in LAL frame */
  REAL8 *S1y,             /**< [out] Dimensionless spin1y in LAL frame */
  REAL8 *S1z,             /**< [out] Dimensionless spin1z in LAL frame */
  REAL8 *S2x,             /**< [out] Dimensionless spin2x in LAL frame */
  REAL8 *S2y,             /**< [out] Dimensionless spin2y in LAL frame */
  REAL8 *S2z,             /**< [out] Dimensionless spin2z in LAL frame */
  REAL8 fRef,             /**< Reference frequency */
  REAL8 mTot,             /**< Total mass */
  const char *NRDataFile  /**< Location of NR HDF file */
);

int XLALSimInspiralNRWaveformGetHplusHcross(
        REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 inclination,              /**< inclination angle */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 r,                        /**< distance of source (m) */
        REAL8 fStart,                   /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z,                      /**< initial value of S2z */
        const char *NRDataFile,         /**< Location of NR HDF file */
        LALValue* ModeArray             /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */
        );

/* in module LALSimIMRNRSur7dq2.c */

double XLALSimInspiralNRSur7dq2StartFrequency(
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z                      /**< initial value of S2z */
);

int XLALSimInspiralNRSur7dq2Polarizations(
        REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 inclination,              /**< inclination angle */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 distnace,                 /**< distance of source (m) */
        REAL8 fMin,                     /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 s1x,                      /**< reference value of S1x */
        REAL8 s1y,                      /**< reference value of S1y */
        REAL8 s1z,                      /**< reference value of S1z */
        REAL8 s2x,                      /**< reference value of S2x */
        REAL8 s2y,                      /**< reference value of S2y */
        REAL8 s2z,                      /**< reference value of S2z */
        LALValue* ModeArray             /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */
);

SphHarmTimeSeries *XLALSimInspiralNRSur7dq2Modes(
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 S1x,                      /**< x-component of the dimensionless spin of object 1 */
        REAL8 S1y,                      /**< y-component of the dimensionless spin of object 1 */
        REAL8 S1z,                      /**< z-component of the dimensionless spin of object 1 */
        REAL8 S2x,                      /**< x-component of the dimensionless spin of object 2 */
        REAL8 S2y,                      /**< y-component of the dimensionless spin of object 2 */
        REAL8 S2z,                      /**< z-component of the dimensionless spin of object 2 */
        REAL8 fMin,                     /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 distnace,                 /**< distance of source (m) */
        LALValue* ModeArray             /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */
);

/* in module LALSimNRTunedTides.c */
double XLALSimNRTunedTidesComputeKappa2T(
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
);

double XLALSimNRTunedTidesMergerFrequency(
    const REAL8 mtot_MSUN, /**< total mass of system (solar masses) */
    const REAL8 kappa2T,   /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
    const REAL8 q          /**< mass-ratio q >= 1 */
);

int XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
    const REAL8Sequence *phi_tidal, /**< [out] tidal phase frequency series */
    const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    );

/* In module LALSimIMRPhenomD_NRTidal.c */
int XLALSimIMRPhenomDNRTidal(COMPLEX16FrequencySeries **htilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, LALDict *extraParams);
int XLALSimIMRPhenomDNRTidalFrequencySequence(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, LALDict *extraParams);

/* in module LALSimIMRPhenomHM.c */
int XLALSimIMRPhenomHM(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 chi1z,
    REAL8 chi2z,
    const REAL8 distance,
    const REAL8 inclination,
    const REAL8 phiRef,
    const REAL8 deltaF,
    REAL8 f_ref,
    LALDict *extraParams);

int XLALSimIMRPhenomHMGethlmModes(
    SphHarmFrequencySeries **hlms,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 chi1z,
    REAL8 chi2z,
    const REAL8 phiRef,
    const REAL8 deltaF,
    REAL8 f_ref,
    LALDict *extraParams);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMR_H */
