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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
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
 * @defgroup LALSimIMRPhenomX_c                  LALSimIMRPhenomX.c
 * @defgroup LALSimIMRPhenomT_c                  LALSimIMRPhenomTHM.c
 * @defgroup LALSimIMREOBNRv2_c                  LALSimIMREOBNRv2.c
 * @defgroup LALSimIMRSpinAlignedEOB_c           LALSimIMRSpinAlignedEOB.c
 * @defgroup LALSimIMRSpinPrecEOB_c              LALSimIMRSpinPrecEOB.c
 * @defgroup LALSimIMRSpinPrecEOBv4P_c           LALSimIMRSpinPrecEOBv4P.c
 * @defgroup LALSimIMRSEOBNRROM_c                LALSimIMRSEOBNRvxROMXXX.c
 * @defgroup LALSimIMRSEOBNRHMROM_c              LALSimIMRSEOBNRv4HMROM.c
 * @defgroup LALSimIMRSEOBNRv2ChirpTime_c        LALSimIMRSEOBNRv2ChirpTime.c
 * @defgroup LALSimIMRPSpinInspiralRD_c          LALSimIMRPSpinInspiralRD.c
 * @defgroup LALSimIMRTidal_c                    LALSimIMRLackeyTidal2013.c
 * @defgroup LALSimPrecessingNRSur_c             LALSimIMRPrecessingNRSur.c
 * @defgroup LALSimIMRNRWaveforms_c              LALSimIMRNRWaveforms.c
 * @defgroup LALSimIMRTEOBResumS_c               LALSimIMRTEOBResumS.c
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
 IMRPhenomPv2NRTidal_V, /**< version Pv2_NRTidal: based on IMRPhenomPv2; NRTides added before precession; can be used with both NRTidal versions defined below */
 IMRPhenomPv3_V  /**< version 3: based on IMRPhenomD and the precession angles from Katerina Chatziioannou PhysRevD.95.104004 (arxiv:1703.03967) */
} IMRPhenomP_version_type;

typedef enum tagNRTidal_version_type {
 NRTidal_V, /**< version NRTidal: based on https://arxiv.org/pdf/1706.02969.pdf*/
 NRTidalv2_V, /**< version NRTidalv2: https://arxiv.org/abs/1905.06011 */
 NRTidalv2NoAmpCorr_V, /**< version NRTidalv2, without amplitude corrections */
 NRTidalv2NSBH_V, /**< version NRTidalv2: https://arxiv.org/abs/1905.06011 with amplitude corrections for NSBH (used for SEOBNRv4ROM_NRTidalv2_NSBH) */
 NoNRT_V /**< special case for PhenomPv2 BBH baseline */
} NRTidal_version_type;

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

/* in module LALSimIMRPhenomNSBH.c */
int XLALSimIMRPhenomNSBHProperties(
    REAL8 *f_RD,                        /**< Output: NSBH ringdown frequency [Hz] */
    REAL8 *f_tide,                      /**< Output: NSBH tidal disruption frequency [Hz] */
    REAL8 *torus_mass,                  /**< Output: Torus remnant mass (kg) */
    REAL8 *compactness,                 /**< Output: Compactness of neutron star */
    REAL8 *final_mass,                  /**< Output: final mass after merger (kg) */
    REAL8 *chif,                        /**< Output: final dimensionless spin */
    REAL8 mBH_SI,                       /**< Mass of BH (kg) */
    REAL8 mNS_SI,                       /**< Mass of neutron star 2 (kg) */
    REAL8 chi_BH,                        /**< Dimensionless aligned component spin of Black Hole */
    REAL8 lambda_NS                     /**< Dimensionless tidal deformability of NS */
);

double XLALSimIMRPhenomNSBH_x_D(const REAL8 Mtorus, const REAL8 C, const REAL8 q, const REAL8 chi);
double XLALSimIMRPhenomNSBH_epsilon_ins_with_torus_mass(const REAL8 Mtorus, const REAL8 C, const REAL8 q, const REAL8 chi);
double XLALSimIMRPhenomNSBH_x_D_prime(const REAL8 Mtorus, const REAL8 C, const REAL8 q, const REAL8 chi);
double XLALSimIMRPhenomNSBH_sigma_tide_with_torus_mass(const REAL8 Mtorus, const REAL8 C, const REAL8 q, const REAL8 chi);
double XLALSimIMRPhenomNSBH_epsilon_tide_ND(const REAL8 x_ND);
double XLALSimIMRPhenomNSBH_sigma_tide_ND(const REAL8 x_ND_prime);
double XLALSimIMRPhenomNSBH_x_ND(const REAL8 f_tide, const REAL8 f_RD_tilde, const REAL8 C, const REAL8 chi);
double XLALSimIMRPhenomNSBH_x_ND_prime(const REAL8 f_tide, const REAL8 f_RD_tilde, const REAL8 C, const REAL8 chi);
double XLALSimIMRPhenomNSBH_delta2_prime(const REAL8 f_tide, const REAL8 f_RD_tilde);

double XLALSimIMRPhenomNSBH_window_plus(const REAL8 f, const REAL8 f0, const REAL8 d);
double XLALSimIMRPhenomNSBH_window_minus(const REAL8 f, const REAL8 f0, const REAL8 d);
double XLALSimIMRPhenomNSBH_eta_from_q(const REAL8 q);

double XLALSimIMRPhenomNSBH_baryonic_mass_from_C(const REAL8 C, const REAL8 Mg);

COMPLEX16 XLALSimIMRPhenomNSBH_omega_tilde(const REAL8 a);

/* in module LALSimNSBHProperties.c */
double XLALSimNSBH_fGWinKerr(const REAL8 r, const REAL8 M, const REAL8 a);
double XLALSimNSBH_rKerrISCO(const REAL8 a);
double XLALSimNSBH_xi_tide(const REAL8 q, const REAL8 a, const REAL8 mu);
double XLALSimNSBH_compactness_from_lambda(const REAL8 Lambda);
double XLALSimNSBH_torus_mass_fit(const REAL8 q, const REAL8 a, const REAL8 C);

/* in module LALSimIMRPhenomD.c */
int XLALSimIMRPhenomDGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 fRef, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi1, const REAL8 chi2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, LALDict *extraParams, NRTidal_version_type NRTidal_version);
int XLALSimIMRPhenomDFrequencySequence(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi0, const REAL8 fRef_in, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi1, const REAL8 chi2, const REAL8 distance, LALDict *extraParams, NRTidal_version_type NRTidal_version);
double XLALIMRPhenomDGetPeakFreq(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in);
double XLALSimIMRPhenomDChirpTime(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in, const REAL8 fHz);
double XLALSimIMRPhenomDFinalSpin(const REAL8 m1_in, const REAL8 m2_in, const REAL8 chi1_in, const REAL8 chi2_in);

int XLALSimIMRPhenomP(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 chi1_l, const REAL8 chi2_l, const REAL8 chip, const REAL8 thetaJ, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, const REAL8 f_ref, IMRPhenomP_version_type IMRPhenomP_version, NRTidal_version_type NRTidal_version, LALDict *extraParams);
int XLALSimIMRPhenomPFrequencySequence(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, const REAL8 chi1_l, const REAL8 chi2_l, const REAL8 chip, const REAL8 thetaJ, REAL8 m1_SI, const REAL8 m2_SI, const REAL8 distance, const REAL8 alpha0, const REAL8 phic, const REAL8 f_ref, IMRPhenomP_version_type IMRPhenomP_version, NRTidal_version_type NRTidal_version, LALDict *extraParams);
int XLALSimIMRPhenomPCalculateModelParametersOld(REAL8 *chi1_l, REAL8 *chi2_l, REAL8 *chip, REAL8 *thetaJ, REAL8 *alpha0, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 lnhatx, const REAL8 lnhaty, const REAL8 lnhatz, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, IMRPhenomP_version_type IMRPhenomP_version);
int XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(REAL8 *chi1_l, REAL8 *chi2_l, REAL8 *chip, REAL8 *thetaJN, REAL8 *alpha0, REAL8 *phi_aligned, REAL8 *zeta_polariz, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 phiRef, const REAL8 incl, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, IMRPhenomP_version_type IMRPhenomP_version);

/* in module LALSimIMREOBNRv2.c */

int XLALSimIMREOBNRv2DominantMode(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
int XLALSimIMREOBNRv2AllModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 phiC, const REAL8 deltaT, const REAL8 m1SI, const REAL8 m2SI, const REAL8 fLower, const REAL8 distance, const REAL8 inclination);
SphHarmTimeSeries *XLALSimIMREOBNRv2Modes(const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 fLower, const REAL8 distance);


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

int XLALSimIMRSEOBNRv4ROM(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, LALDict *LALparams, NRTidal_version_type NRTidal_version);
int XLALSimIMRSEOBNRv4ROMFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, LALDict *LALparams, NRTidal_version_type NRTidal_version);
int XLALSimIMRSEOBNRv4ROMTimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
int XLALSimIMRSEOBNRv4ROMFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);

/* in module LALSimIMRSEOBNRv4HMROM.c */
int XLALSimIMRSEOBNRv4HMROM(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, UINT4 nModes, bool use_hybridization, LALDict *LALParams);
int XLALSimIMRSEOBNRv4HMROMFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, UINT4 nModes, LALDict *LALParams);
int XLALSimIMRSEOBNRv4HMROM_Modes(SphHarmFrequencySeries **hlm, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, UINT4 nModes, bool use_hybridization);
int XLALSimIMRSEOBNRv4HMROMFrequencySequence_Modes(SphHarmFrequencySeries **hlm, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, INT4 nk_max, UINT4 nModes, LALDict *LALParams);


/* in module LALSimIMRSEOBNRv4ROM_NRTidal.c */

int XLALSimIMRSEOBNRv4ROMNRTidalFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 Lambda1, REAL8 Lambda2, LALDict *LALparams, NRTidal_version_type NRTidal_version);
int XLALSimIMRSEOBNRv4ROMNRTidal(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 Lambda1, REAL8 Lambda2, LALDict *LALparams, NRTidal_version_type NRTidal_version);

/* in module LALSimBHNSRemnantFits.c */
REAL8 XLALbbh_final_mass_non_precessing_UIB2016(const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 chi2);
REAL8 XLALbbh_final_spin_non_precessing_UIB2016(const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 chi2);
REAL8 XLALBHNS_mass_aligned(const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 lam);
REAL8 XLALBHNS_spin_aligned(const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 lam);

/* In module LALSimIMRSEOBNRv4ROM_NSBHAmplitudeCorrection.c */
int XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries(
    const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 chi1, /**< Spin of black hole */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
);

/* in module LALSimIMRSEOBNRv4TSurrogate.c */

int XLALSimIMRSEOBNRv4TSurrogate(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, SEOBNRv4TSurrogate_spline_order spline_order);
int XLALSimIMRSEOBNRv4TSurrogateFrequencySequence(struct tagCOMPLEX16FrequencySeries **hptilde, struct tagCOMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, SEOBNRv4TSurrogate_spline_order spline_order);


/* in module LALSimIMRPSpinInspiralRD.c */

int XLALSimIMRPhenSpinFinalMassSpin(REAL8 *finalMass, REAL8 *finalSpin, REAL8 m1, REAL8 m2, REAL8 s1s1, REAL8 s2s2, REAL8 s1L, REAL8 s2L, REAL8 s1s2, REAL8 energy);
int XLALSimSpinInspiralGenerator(REAL8TimeSeries **hPlus, REAL8TimeSeries **hCross, REAL8 phi_start, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams);
int XLALSimIMRPhenSpinInspiralRDGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phi0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 iota, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, int phaseO, int ampO, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams);

/* IMRPhenomX/HM Routines */
/* in module LALSimIMRPhenomX.c */
int XLALSimIMRPhenomXASGenerateFD(COMPLEX16FrequencySeries **htilde22,
  REAL8 m1_SI,
  REAL8 m2_SI,
  REAL8 chi1L,
  REAL8 chi2L,
  REAL8 distance,
  REAL8 f_min,
  REAL8 f_max,
  REAL8 deltaF,
  REAL8 phiRef,
  REAL8 fRef_In,
  LALDict *lalParams
);

int XLALSimIMRPhenomXASFrequencySequence(
  COMPLEX16FrequencySeries **htilde22,
  const REAL8Sequence *freqs,
  REAL8 m1_SI,
  REAL8 m2_SI,
  REAL8 chi1L,
  REAL8 chi2L,
  REAL8 distance,
  REAL8 phiRef,
  REAL8 fRef_In,
  LALDict *lalParams
);

int XLALSimIMRPhenomXPMSAAngles(
 REAL8Sequence *phiz_of_f,                 /**< [out]  */
 REAL8Sequence *zeta_of_f,                 /**< [out] */
 REAL8Sequence *costhetaL_of_f,            /**< [out] */
 const REAL8Sequence *freqs,               /**< [out] Frequency series [Hz]         */
 const REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
 const REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
 const REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 fRef_In,                      /**< Reference frequency (Hz) */
 LALDict *lalParams                        /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPPNAngles(
 REAL8Sequence *alpha_of_f,          /**< [out]  */
 REAL8Sequence *gamma_of_f,          /**< [out] */
 REAL8Sequence *cosbeta_of_f,        /**< [out] */
 const REAL8Sequence *freqs,         /**< [out] Frequency series [Hz]         */
 const REAL8 m1_SI,                  /**< mass of companion 1 (kg) */
 const REAL8 m2_SI,                  /**< mass of companion 2 (kg) */
 const REAL8 chi1x,                  /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi1y,                  /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi1z,                  /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2x,                  /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2y,                  /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 chi2z,                  /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
 const REAL8 fRef_In,                /**< Reference frequency (Hz) */
 LALDict *lalParams                  /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde,         /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,         /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
  REAL8 chi1x,                                /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                                /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                                /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                                /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                                /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                                /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  const REAL8 distance,                       /**< Distance of source (m) */
  const REAL8 inclination,                    /**< inclination of source (rad) */
  const REAL8 phiRef,                         /**< Orbital phase (rad) at reference frequency */
  REAL8 f_min,                                /**< Starting GW frequency (Hz) */
  REAL8 f_max,                                /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
  const REAL8 deltaF,                         /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
  REAL8 fRef_In,                              /**< Reference frequency (Hz) */
  LALDict *lalParams                          /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPFrequencySequence(
  COMPLEX16FrequencySeries **hptilde,         /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,         /**< [out] Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                 /**< [out] Frequency series [Hz]         */
  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
  REAL8 chi1x,                                /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                                /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                                /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                                /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                                /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                                /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  const REAL8 distance,                       /**< Distance of source (m) */
  const REAL8 inclination,                    /**< inclination of source (rad) */
  const REAL8 phiRef,                         /**< Orbital phase (rad) at reference frequency */
  REAL8 fRef_In,                              /**< Reference frequency (Hz) */
  LALDict *lalParams                          /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXHMModes(
     SphHarmFrequencySeries **hlms,
  	  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
      REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
      REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
      REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
      REAL8 deltaF,                               /**< sampling interval (s) */
  		REAL8 f_min,                                /**< starting GW frequency (Hz) */
  		REAL8 f_max,                                /**< ending GW frequency (Hz) */
      REAL8 f_ref,                                /**< reference GW frequency (Hz) */
      REAL8 phiRef,                               /**< phase shift at reference frequency */
      REAL8 distance,                             /**< distance of source (m) */
      LALDict *LALparams                          /**< LAL dictionary with extra options */
);

int XLALSimIMRPhenomXHMGenerateFDOneMode(
 COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
 REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
 REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
 REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
 REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
 UINT4 ell,                           /**< l index of the mode */
 INT4 emm,                            /**< m index of the mode */
 REAL8 distance,                      /**< Luminosity distance (m) */
 REAL8 f_min,                         /**< Starting GW frequency (Hz) */
 REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
 REAL8 deltaF,                        /**< Sampling frequency (Hz) */
 REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
 REAL8 fRef_In,                       /**< Reference frequency (Hz) */
 LALDict *lalParams                   /**< lal dictionary parameters */
);

int XLALSimIMRPhenomXHMFrequencySequenceOneMode(
    COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
    REAL8Sequence *freqs,                /**< frequency array to evaluate model */
    REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
    REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
    REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
    REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
    UINT4 ell,                           /**< l index of the mode */
    INT4 emm,                            /**< m index of the mode */
    REAL8 distance,                      /**< Luminosity distance (m) */
    REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
    REAL8 fRef_In,                       /**< Reference frequency (Hz) */
    LALDict *lalParams                   /**< lal dictionary parameters */
);

int XLALSimIMRPhenomXHM(
   COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
   COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
   REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
   REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
   REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
   REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
   REAL8 f_min,                        /**< Starting GW frequency (Hz) */
   REAL8 f_max,                        /**< End frequency; 0 defaults to Mf = 0.3 */
   REAL8 deltaF,                       /**< Sampling frequency (Hz) */
   REAL8 distance,                     /**< distance of source (m) */
   REAL8 inclination,                  /**< inclination of source (rad) */
   REAL8 phiRef,                       /**< reference orbital phase (rad) */
   REAL8 fRef_In,                      /**< Reference frequency */
   LALDict *lalParams                  /**<linked list containing the extra testing GR parameters */
);


int XLALSimIMRPhenomXHM2(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 inclination,                   /**  Inclination of the source */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< LAL Dictionary */
);


int XLALSimIMRPhenomXHMMultiBandOneMode(
  COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  UINT4 ell,                           /**< l index of the mode */
  INT4 emm,                            /**< m index of the mode */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< Extra params */
);

int XLALSimIMRPhenomXHMMultiBandOneModeMixing(
  COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
  COMPLEX16FrequencySeries *htilde22,  /**< [out] FD waveform */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  UINT4 ell,                           /**< l index of the mode */
  INT4 emm,                            /**< m index of the mode */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< Extra params */
);

int XLALSimIMRPhenomXHMAmplitude(
    REAL8FrequencySeries **amplitude, /**< [out] FD amp */
    REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
    REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
    REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
    REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
    UINT4 ell,                           /**< l index of the mode */
    INT4 emm,                            /**< m index of the mode */
    REAL8 distance,                      /**< Luminosity distance (m) */
    REAL8 f_min,                         /**< Starting GW frequency (Hz) */
    REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
    REAL8 deltaF,                        /**< Sampling frequency (Hz) */
    REAL8 phiRef,                        /**< Orbital amp at fRef (rad) */
    REAL8 fRef_In,                       /**< Reference frequency (Hz) */
    LALDict *lalParams                   /**< Extra params */
  );

int XLALSimIMRPhenomXHMPhase(
    REAL8FrequencySeries **phase,        /**< [out] FD amp */
    REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
    REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
    REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
    REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
    UINT4 ell,                           /**< l index of the mode */
    INT4 emm,                            /**< m index of the mode */
    REAL8 distance,                      /**< Luminosity distance (m) */
    REAL8 f_min,                         /**< Starting GW frequency (Hz) */
    REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
    REAL8 deltaF,                        /**< Sampling frequency (Hz) */
    REAL8 phiRef,                        /**< Orbital amp at fRef (rad) */
    REAL8 fRef_In,                       /**< Reference frequency (Hz) */
    LALDict *lalParams                   /**< Extra params */
  );

int XLALSimIMRPhenomXPHM(
  COMPLEX16FrequencySeries **hptilde,         /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,         /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
  REAL8 chi1x,                                /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                                /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                                /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                                /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                                /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                                /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 distance,                             /**< Distance of source (m) */
  REAL8 inclination,                          /**< inclination of source (rad) */
  REAL8 phiRef,                               /**< Orbital phase (rad) at reference frequency */
  REAL8 f_min,                                /**< Starting GW frequency (Hz) */
  REAL8 f_max,                                /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
  REAL8 deltaF,                               /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
  REAL8 fRef_In,                              /**< Reference frequency (Hz) */
  LALDict *lalParams                          /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPHMFromModes(
COMPLEX16FrequencySeries **hptilde,         /**< [out] Frequency-domain waveform h+ */
COMPLEX16FrequencySeries **hctilde,         /**< [out] Frequency-domain waveform hx */
REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
REAL8 chi1x,                                /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
REAL8 chi1y,                                /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
REAL8 chi1z,                                /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
REAL8 chi2x,                                /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
REAL8 chi2y,                                /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
REAL8 chi2z,                                /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
const REAL8 distance,                       /**< Distance of source (m) */
const REAL8 inclination,                    /**< inclination of source (rad) */
const REAL8 phiRef,                         /**< Orbital phase (rad) at reference frequency */
REAL8 f_min,                                /**< Starting GW frequency (Hz) */
REAL8 f_max,                                /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
const REAL8 deltaF,                         /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
REAL8 fRef_In,                              /**< Reference frequency (Hz) */
LALDict *lalParams                          /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXHMFrequencySequence(
   COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+  */
   COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx  */
   REAL8Sequence *freqs,               /**< Input Frequency series [Hz]         */
   REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
   REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
   REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
   REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
   REAL8 distance,                     /**< Distance of source (m) */
   REAL8 inclination,                  /**< inclination of source (rad) */
   REAL8 phiRef,                       /**< Orbital phase (rad) at reference frequency */
   REAL8 fRef,                         /**< Reference frequency (Hz) */
   LALDict *lalParams                  /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPHMFrequencySequence(
  COMPLEX16FrequencySeries **hptilde,         /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,         /**< [out] Frequency-domain waveform hx */
  REAL8Sequence *freqs,                       /**< input frequency series [Hz]         */
  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
  REAL8 chi1x,                                /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                                /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                                /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                                /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                                /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                                /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 distance,                             /**< Distance of source (m) */
  REAL8 inclination,                          /**< inclination of source (rad) */
  REAL8 phiRef,                               /**< Orbital phase (rad) at reference frequency */
  REAL8 fRef_In,                              /**< Reference frequency (Hz) */
  LALDict *lalParams                          /**< LAL Dictionary struct */
);

int XLALSimIMRPhenomXPHMOneMode(
  COMPLEX16FrequencySeries **hlmpos,      /**< [out] Frequency-domain waveform hlm inertial frame positive frequencies */
  COMPLEX16FrequencySeries **hlmneg,      /**< [out] Frequency-domain waveform hlm inertial frame negative frequencies */
  const UINT4 l,                          /**< First index of the (l,m) precessing mode */
  const INT4  m,                          /**< Second index of the (l,m) precessing mode */
  const REAL8 m1_SI,                      /**< mass of companion 1 (kg) */
  const REAL8 m2_SI,                      /**< mass of companion 2 (kg) */
  const REAL8 chi1x,                      /**< x-component of the dimensionless spin of object 1 */
  const REAL8 chi1y,                      /**< y-component of the dimensionless spin of object 1 */
  const REAL8 chi1z,                      /**< z-component of the dimensionless spin of object 1 */
  const REAL8 chi2x,                      /**< x-component of the dimensionless spin of object 2 */
  const REAL8 chi2y,                      /**< y-component of the dimensionless spin of object 2 */
  const REAL8 chi2z,                      /**< z-component of the dimensionless spin of object 2 */
  const REAL8 distance,                   /**< distance of source (m) */
  const REAL8 inclination,                /**< inclination of source (rad) */
  const REAL8 phiRef,                     /**< reference orbital phase (rad) */
  const REAL8 deltaF,                     /**< Sampling frequency (Hz) */
  const REAL8 f_min,                      /**< Starting GW frequency (Hz) */
  const REAL8 f_max,                      /**< End frequency; 0 defaults to ringdown cutoff freq */
  const REAL8 fRef_In,                    /**< Reference frequency */
  LALDict *lalParams                      /**<LAL Dictionary */
);

int XLALSimIMRPhenomXPHMModes(
      SphHarmFrequencySeries **hlms,              /**< [out] list with single modes h_lm in the J-frame */
  	  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
      REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
      REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 deltaF,                               /**< frequency spacing (Hz) */
  		REAL8 f_min,                                /**< starting GW frequency (Hz) */
  		REAL8 f_max,                                /**< ending GW frequency (Hz) */
      REAL8 f_ref,                                /**< reference GW frequency (Hz) */
      REAL8 phiRef,                               /**< phase shift at reference frequency */
      REAL8 distance,                             /**< distance of source (m) */
  	  REAL8 inclination,                          /**< inclination of source (rad) */
    	LALDict *LALparams                          /**< LAL dictionary with extra options */
);

/*
		XLAL routine to calculate the IMRPhenomX mode parameters from the source frame.
		Note that this is just an IMRPhenomX compatible implementation of:

		XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame
*/
int XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame(
    REAL8 *chi1L,                     /**< [out] Dimensionless aligned spin on companion 1 */
    REAL8 *chi2L,                     /**< [out] Dimensionless aligned spin on companion 2 */
    REAL8 *chi_p,                     /**< [out] Effective spin in the orbital plane */
    REAL8 *thetaJN,                   /**< [out] Angle between J0 and line of sight (z-direction) */
    REAL8 *alpha0,                    /**< [out] Initial value of alpha angle (azimuthal precession angle) */
    REAL8 *phi_aligned,               /**< [out] Initial phase to feed the underlying aligned-spin model */
    REAL8 *zeta_polarization,         /**< [out] Angle to rotate the polarizations */
    const REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
    const REAL8 f_ref,                /**< Reference GW frequency (Hz) */
    const REAL8 phiRef,               /**< Reference phase */
    const REAL8 incl,                 /**< Inclination : angle between LN and the line of sight */
    const REAL8 s1x,                  /**< Initial value of s1x: dimensionless spin of BH 1 */
    const REAL8 s1y,                  /**< Initial value of s1y: dimensionless spin of BH 1 */
    const REAL8 s1z,                  /**< Initial value of s1z: dimensionless spin of BH 1 */
    const REAL8 s2x,                  /**< Initial value of s2x: dimensionless spin of BH 2 */
    const REAL8 s2y,                  /**< Initial value of s2y: dimensionless spin of BH 2 */
    const REAL8 s2z,                  /**< Initial value of s2z: dimensionless spin of BH 2 */
    LALDict *lalParams                /**< LAL Dictionary */
);

/* IMRPhenomT/HM Routines */
/* in module LALSimIMRPhenomTHM.c */

SphHarmTimeSeries *XLALSimIMRPhenomTHM_Modes(
  REAL8 m1_SI,        /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,        /**< Mass of companion 2 (kg) */
  REAL8 chi1L,        /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,        /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,     /**< Luminosity distance (m) */
  REAL8 deltaT,       /**< inclination of source (rad) */
  REAL8 fmin,         /**< sampling interval (s) */
  REAL8 fRef,         /**< reference GW frequency (Hz) */
  REAL8 phiRef,       /**< reference orbital phase (rad) */
  LALDict *lalParams /**< LAL dictionary containing accessory parameters */
  );

int XLALSimIMRPhenomTHM(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,          /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,          /**< Mass of companion 2 (kg) */
  REAL8 chi1L,          /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,          /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,       /**< Luminosity distance (m) */
  REAL8 inclination,    /**< inclination of source (rad) */
  REAL8 deltaT,         /**< sampling interval (s) */
  REAL8 fmin,           /**< starting GW frequency (Hz) */
  REAL8 fRef,           /**< reference GW frequency (Hz) */
  REAL8 phiRef,         /**< reference orbital phase (rad) */
  LALDict *lalparams  /**< LAL dictionary containing accessory parameters */
  );

int XLALSimIMRPhenomT(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1L,      /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,      /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  );

/* IMRPhenomTPHM Routines */
/* in module LALSimIMRPhenomTPHM.c */

int XLALSimIMRPhenomTP(
  REAL8TimeSeries **hp,       /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc,       /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams            /**< LAL dictionary containing accessory parameters */
  );

int XLALSimIMRPhenomTPHM(
  REAL8TimeSeries **hp,       /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc,       /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams            /**< LAL dictionary containing accessory parameters */
  );

SphHarmTimeSeries* XLALSimIMRPhenomTPHM_ChooseTDModes(
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  LALDict *lalParams            /**< LAL dictionary containing accessory parameters */
  );

int XLALSimIMRPhenomTPHM_L0Modes(
  SphHarmTimeSeries **hlmI,   /**< [out] Modes in the intertial L0=z frame*/
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,            /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  );

int XLALSimIMRPhenomTPHM_JModes(
  SphHarmTimeSeries **hlmJ,   /**< [out] Modes in the intertial J0=z frame*/
  REAL8TimeSeries **alphaTS,  /**< [out] Precessing Euler angle alpha */
  REAL8TimeSeries **cosbetaTS,   /**< [out] cosinus of Precessing Euler angle beta */
  REAL8TimeSeries **gammaTS,  /**< [out] Precessing Euler angle gamma */
  REAL8 *af,                  /**< [out] Final spin */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,       /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  );

int XLALSimIMRPhenomTPHM_EvolveOrbit(
  REAL8TimeSeries **V,            /**< post-Newtonian parameter [returned]*/
  REAL8TimeSeries **S1x,          /**< Spin1 vector x component [returned]*/
  REAL8TimeSeries **S1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S1z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **S2x,          /**< Spin2 vector x component [returned]*/
  REAL8TimeSeries **S2y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S2z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **LNhatx,       /**< unit orbital ang. mom. x [returned]*/
  REAL8TimeSeries **LNhaty,       /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **LNhatz,       /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **E1x,          /**< orb. plane basis vector x[returned]*/
  REAL8TimeSeries **E1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **E1z,          /**< "    "    "  z component [returned]*/
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams       /**< LAL dictionary containing accessory parameters */
  );

int XLALSimIMRPhenomTPHM_CoprecModes(
  SphHarmTimeSeries **hlmJ,   /**< [out] Modes in the intertial J0=z frame*/
  REAL8TimeSeries **alphaTS,  /**< [out] Precessing Euler angle alpha */
  REAL8TimeSeries **cosbetaTS,   /**< [out] Precessing Euler angle beta */
  REAL8TimeSeries **gammaTS,  /**< [out] Precessing Euler angle gamma */
  REAL8 *af,                  /**< [out] Final spin */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,       /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  );

/* in module LALSimIMRTEOBResumS.c */

int XLALSimIMRTEOBResumS(
                         REAL8TimeSeries **hplus,
                         REAL8TimeSeries **hcross,
                         const REAL8 phiRef,
                         const REAL8 deltaT,
                         const REAL8 m1,
                         const REAL8 m2,
                         const REAL8 S1x,
                         const REAL8 S1y,
                         const REAL8 S1z,
                         const REAL8 S2x,
                         const REAL8 S2y,
                         const REAL8 S2z,
                         const REAL8 lambda1,
                         const REAL8 lambda2,
                         const REAL8 distance,
                         const REAL8 inclination,
                         const REAL8 longAscNodes,
                         LALDict *LALparams,
                         const REAL8 eccentricity,
                         const REAL8 meanPerAno,
                         const REAL8 f_min,
                         const REAL8 f_ref
                         );


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

/* The following XLALSimInspiralNRWaveformGetHplusHcross() generates polarizations
 * reading directly the NR files and does not return l,m modes.
 */
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

/* The following XLALSimInspiralNRWaveformGetHlms() reads NR file to output l,m modes.
 */
INT4 XLALSimInspiralNRWaveformGetHlms(SphHarmTimeSeries **hlms, /**< OUTPUT */
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

/* in module LALSimIMRPrecessingNRSur.c */

int XLALSimInspiralPrecessingNRSurPolarizations(
        REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< azimuthal angle for Ylms */
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
        LALDict* LALparams,             /**< Dict with extra parameters */
        Approximant approximant     /**< approximant (NRSur7dq2 or NRSur7dq4) */

);

SphHarmTimeSeries *XLALSimInspiralPrecessingNRSurModes(
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
        REAL8 distance,                 /**< distance of source (m) */
        LALDict* LALparams,             /**< Dict with extra parameters */
        Approximant approximant     /**< approximant (NRSur7dq2 or NRSur7dq4) */
);

int XLALPrecessingNRSurDynamics(
        gsl_vector **t_dynamics, /**< Output: Time array at which the dynamics are returned. */
        gsl_vector **quat0,      /**< Output: Time series of 0th index of coprecessing frame quaternion. */
        gsl_vector **quat1,      /**< Output: Time series of 1st index of coprecessing frame quaternion. */
        gsl_vector **quat2,      /**< Output: Time series of 2nd index of coprecessing frame quaternion. */
        gsl_vector **quat3,      /**< Output: Time series of 3rd index of coprecessing frame quaternion. */
        gsl_vector **orbphase,   /**< Output: Time series of orbital phase in the coprecessing frame. */
        gsl_vector **chiAx,      /**< Output: Time series of x-comp of dimensionless spin of BhA in the coprecessing frame. */
        gsl_vector **chiAy,      /**< Output: Time series of y-comp of dimensionless spin of BhA in the coprecessing frame. */
        gsl_vector **chiAz,      /**< Output: Time series of z-comp of dimensionless spin of BhA in the coprecessing frame. */
        gsl_vector **chiBx,      /**< Output: Time series of x-comp of dimensionless spin of BhB in the coprecessing frame. */
        gsl_vector **chiBy,      /**< Output: Time series of y-comp of dimensionless spin of BhB in the coprecessing frame. */
        gsl_vector **chiBz,      /**< Output: Time series of z-comp of dimensionless spin of BhB in the coprecessing frame. */
        REAL8 q,                 /**< mass ratio m1/m2 >= 1. */
        REAL8 chiA0x,            /**< x-comp of dimensionless spin of BhA in the coorbital frame at the reference epoch. */
        REAL8 chiA0y,            /**< y-comp of dimensionless spin of BhA in the coorbital frame at the reference epoch. */
        REAL8 chiA0z,            /**< z-comp of dimensionless spin of BhA in the coorbital frame at the reference epoch. */
        REAL8 chiB0x,            /**< x-comp of dimensionless spin of BhB in the coorbital frame at the reference epoch. */
        REAL8 chiB0y,            /**< y-comp of dimensionless spin of BhB in the coorbital frame at the reference epoch. */
        REAL8 chiB0z,            /**< z-comp of dimensionless spin of BhB in the coorbital frame at the reference epoch. */
        REAL8 omegaRef_dimless,  /**< Dimensionless orbital frequency (rad/M) in the coprecessing frame at the reference epoch.*/
        REAL8 init_quat0,        /**< 0th comp of the coprecessing frame quaternion at the reference epoch.*/
        REAL8 init_quat1,        /**< 1st comp of the coprecessing frame quaternion at the reference epoch.*/
        REAL8 init_quat2,        /**< 2nd comp of the coprecessing frame quaternion at the reference epoch.*/
        REAL8 init_quat3,        /**< 3rd comp of the coprecessing frame quaternion at the reference epoch.*/
        REAL8 init_orbphase,     /**< orbital phase in the coprecessing frame at the reference epoch. */
        LALDict* LALparams,      /**< Dict with extra parameters. */
        Approximant approximant  /**< approximant (NRSur7dq2 or NRSur7dq4). */
);

/* in module LALSimNRSur7dq4Remnant.c */
int XLALNRSur7dq4Remnant(
    gsl_vector **result,        /**<Output: The requested remnant property. */
    REAL8 q,                    /**< Mass ratio of Bh1/Bh2. q>=1. */
    REAL8 s1x,                  /**< S1x in coorbital frame at t=-100M */
    REAL8 s1y,                  /**< S1y in coorbital frame at t=-100M */
    REAL8 s1z,                  /**< S1z in coorbital frame at t=-100M */
    REAL8 s2x,                  /**< S2x in coorbital frame at t=-100M */
    REAL8 s2y,                  /**< S2y in coorbital frame at t=-100M */
    REAL8 s2z,                  /**< S2z in coorbital frame at t=-100M */
    char *remnant_property,     /**< One of "mf", "chif" or "vf" */
    LALDict* LALparams          /**< Dict with extra parameters */
);

/* in module LALSimNRSur3dq8Remnant.c */
int XLALNRSur3dq8Remnant(
    REAL8 *result,              /**<Output: The requested remnant property. */
    REAL8 q,                    /**< Mass ratio of Bh1/Bh2. q>=1. */
    REAL8 s1z,                  /**< S1z z-spin of Bh1 */
    REAL8 s2z,                  /**< S2z z-spin of Bh2 */
    char *remnant_property,     /**< One of "mf", "chifz", "vfx" or "vfy" */
    LALDict* LALparams          /**< Dict with extra parameters */
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

int XLALSimNRTunedTidesFDTidalAmplitudeFrequencySeries(
    const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
);

int XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
    const REAL8Sequence *phi_tidal, /**< [out] tidal phase frequency series */
    const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
    const REAL8Sequence *planck_taper, /**< [out] planck taper */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2, /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    NRTidal_version_type NRTidal_version /**< NRTidal version */
    );

void XLALSimInspiralGetHOSpinTerms(REAL8 *SS_3p5PN, REAL8 *SSS_3p5PN, REAL8 X_A, REAL8 X_B, REAL8 chi1, REAL8 chi2, REAL8 quadparam1, REAL8 quadparam2);

/* In module LALSimIMRPhenomD_NRTidal.c */
int XLALSimIMRPhenomDNRTidal(COMPLEX16FrequencySeries **htilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 fRef, REAL8 distance, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, LALDict *extraParams, NRTidal_version_type NRTidal_version);
int XLALSimIMRPhenomDNRTidalFrequencySequence(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 fRef, REAL8 distance, REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1, REAL8 chi2, REAL8 lambda1, REAL8 lambda2, LALDict *extraParams, NRTidal_version_type NRTidal_version);

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
    REAL8 chi1x,
    REAL8 chi1y,
    REAL8 chi1z,
    REAL8 chi2x,
    REAL8 chi2y,
    REAL8 chi2z,
    const REAL8 phiRef,
    const REAL8 deltaF,
    REAL8 f_ref,
    LALDict *extraParams);

/* from LALSimIMRPhenomNSBH.c */

int XLALSimIMRPhenomNSBHFrequencySequence(
    COMPLEX16FrequencySeries **htilde,
    const REAL8Sequence *freqs,
    REAL8 phiRef,
    REAL8 fRef,
    REAL8 distance,
    REAL8 mBH_SI,
    REAL8 mNS_SI,
    REAL8 chi_BH,
    REAL8 chi_NS,
    LALDict *extraParams);

int XLALSimIMRPhenomNSBH(
    COMPLEX16FrequencySeries **htilde,
    REAL8 phiRef,
    REAL8 deltaF,
    REAL8 fLow,
    REAL8 fHigh,
    REAL8 fRef,
    REAL8 distance,
    REAL8 mBH_SI,
    REAL8 mNS_SI,
    REAL8 chi_BH,
    REAL8 chi_NS,
    LALDict *extraParams
);

/* LALSimInspiralFDPrecAngles functions */
int XLALComputeAngles2PNNonSpinning(
    REAL8Sequence *phiz_of_f,
    REAL8Sequence *zeta_of_f,
    REAL8Sequence *costhetaL_of_f,
    const REAL8Sequence *f,
    const double m1,
    const double m2,
    const double mul,
    const double phl,
    const double mu1,
    const double ph1,
    const double ch1,
    const double mu2,
    const double ph2,
    double ch2,
    const double f_0,
    const int ExpansionOrder);

int XLALComputeAngles3PN(
    REAL8Sequence *phiz_of_f,
    REAL8Sequence *zeta_of_f,
    REAL8Sequence *costhetaL_of_f,
    const REAL8Sequence *f,
    const double m1,
    const double m2,
    const double mul,
    const double phl,
    const double mu1,
    const double ph1,
    const double ch1,
    const double mu2,
    const double ph2,
    double ch2,
    const double f_0,
    const int ExpansionOrder);

int XLALComputeAngles(
    REAL8Sequence *phiz_of_f,
    REAL8Sequence *zeta_of_f,
    REAL8Sequence *costhetaL_of_f,
    const REAL8Sequence *f,
    const double m1,
    const double m2,
    const double mul,
    const double phl,
    const double mu1,
    const double ph1,
    const double ch1,
    const double mu2,
    const double ph2,
    double ch2,
    const double f_0,
    const int ExpansionOrder);

int XLALOrbitalAngMom3PNSpinning(
    REAL8Sequence *L_norm_3PN,
    REAL8Sequence *f_orb_hz,
    const double m1,
    const double m2,
    const double mul,
    const double phl,
    double mu1,
    double ph1,
    double ch1,
    double mu2,
    double ph2,
    double ch2,
    const double f_0,
    const int ExpansionOrder);

/* IMRPhenomPv3 XLAL functions */
int XLALSimIMRPhenomPv3(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 S1x,
    REAL8 S1y,
    REAL8 S1z,
    REAL8 S2x,
    REAL8 S2y,
    REAL8 S2z,
    const REAL8 distance,
    const REAL8 inclination,
    const REAL8 phiRef,
    const REAL8 deltaF,
    const REAL8 f_ref,
    LALDict *extraParams);

/* IMRPhenomPv3HM XLAL functions */
int XLALSimIMRPhenomPv3HMGetHplusHcross(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 chi1x,
    REAL8 chi1y,
    REAL8 chi1z,
    REAL8 chi2x,
    REAL8 chi2y,
    REAL8 chi2z,
    const REAL8 distance,
    const REAL8 inclination,
    const REAL8 phiRef,
    const REAL8 deltaF,
    REAL8 f_ref,
    LALDict *extraParams);

int XLALSimIMRPhenomPv3HMModes(
    SphHarmFrequencySeries **hlms,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 chi1x,
    REAL8 chi1y,
    REAL8 chi1z,
    REAL8 chi2x,
    REAL8 chi2y,
    REAL8 chi2z,
    const REAL8 phiRef,
    const REAL8 deltaF,
    const REAL8 f_ref,
    LALDict *extraParams);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMR_H */
