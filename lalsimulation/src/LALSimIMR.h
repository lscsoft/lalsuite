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
 * @ingroup lalsimulation_general
 *
 * @brief Routines for generating inspiral-merger-ringdown waveforms.
 *
 * @{
 * @defgroup LALSimIMRPhenom_c                LALSimIMRPhenom.c
 * @defgroup LALSimIMREOBNRv2_c               LALSimIMREOBNRv2.c
 * @defgroup LALSimIMRSpinAlignedEOB_c        LALSimIMRSpinAlignedEOB.c
 * @defgroup LALSimIMRSEOBNRv1ROMDoubleSpin_c LALSimIMRSEOBNRv1ROMDoubleSpin.c
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

/*
 * SEOBNRv1 reduced order models
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 */

int XLALSimIMRSEOBNRv1ROMEffectiveSpin(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
    REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
    REAL8 fHigh,                                  /**< End frequency; 0 defaults to ringdown cutoff freq */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi                                     /**< Effective aligned spin */
);

/** Compute waveform in LAL format at specified frequencies */
int XLALSimIMRSEOBNRv1ROMEffectiveSpinFrequencySequence(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi                                     /**< Effective aligned spin */
);

int XLALSimIMRSEOBNRv1ROMDoubleSpin(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
    REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
    REAL8 fHigh,                                  /**< End frequency; 0 defaults to ringdown cutoff freq */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
    REAL8 chi2                                    /**< Dimensionless aligned component spin 2 */
);

/** Compute waveform in LAL format at specified frequencies */
int XLALSimIMRSEOBNRv1ROMDoubleSpinFrequencySequence(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
    REAL8 chi2                                    /**< Dimensionless aligned component spin 2 */
);


/*
 * SEOBNRv2 reduced order models PRELIMINARY!
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 */

int XLALSimIMRSEOBNRv2ROMEffectiveSpin(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
    REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
    REAL8 fHigh,                                  /**< End frequency; 0 defaults to ringdown cutoff freq */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi                                     /**< Effective aligned spin */
);

/** Compute waveform in LAL format at specified frequencies */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi                                     /**< Effective aligned spin */
);

/**
 * Compute the time at a given frequency. The origin of time is at the merger.
 * The allowed frequency range for the input is Mf in [0.0001, 0.3].
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinTimeOfFrequency(
  REAL8 *t,         /**< Output: time (s) at frequency */
  REAL8 frequency,  /**< Frequency (Hz) */
  REAL8 m1SI,       /**< Mass of companion 1 (kg) */
  REAL8 m2SI,       /**< Mass of companion 2 (kg) */
  REAL8 chi         /**< Effective aligned spin */
);

/**
 * Compute the frequency at a given time. The origin of time is at the merger.
 * The frequency range for the output is Mf in [0.0001, 0.3].
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencyOfTime(
  REAL8 *frequency,   /**< Output: Frequency (Hz) */
  REAL8 t,            /**< Time (s) at frequency */
  REAL8 m1SI,         /**< Mass of companion 1 (kg) */
  REAL8 m2SI,         /**< Mass of companion 2 (kg) */
  REAL8 chi           /**< Effective aligned spin */
);

int XLALSimIMRSEOBNRv2ROMDoubleSpin(
    struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
    struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
    REAL8 phiRef,                                 /**< Phase at reference frequency */
    REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
    REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
    REAL8 fHigh,                                  /**< End frequency; 0 defaults to ringdown cutoff freq */
    REAL8 fRef,                                   /**< Reference frequency; 0 defaults to fLow */
    REAL8 distance,                               /**< Distance of source (m) */
    REAL8 inclination,                            /**< Inclination of source (rad) */
    REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
    REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
    REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
    REAL8 chi2                                    /**< Dimensionless aligned component spin 2 */
);

/** Compute waveform in LAL format at specified frequencies */
int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2                                    /**< Dimensionless aligned component spin 2 */
);

 /**
  * Compute the time at a given frequency. The origin of time is at the merger.
  * The allowed frequency range for the input is from Mf = 0.00053 to half the ringdown frequency.
  */
 int XLALSimIMRSEOBNRv2ROMDoubleSpinTimeOfFrequency(
   REAL8 *t,         /**< Output: time (s) at frequency */
   REAL8 frequency,  /**< Frequency (Hz) */
   REAL8 m1SI,       /**< Mass of companion 1 (kg) */
   REAL8 m2SI,       /**< Mass of companion 2 (kg) */
   REAL8 chi1,       /**< Dimensionless aligned component spin 1 */
   REAL8 chi2        /**< Dimensionless aligned component spin 2 */
 );

 /**
  * Compute the frequency at a given time. The origin of time is at the merger.
  * The frequency range for the output is from Mf = 0.00053 to half the ringdown frequency.
  */
 int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencyOfTime(
   REAL8 *frequency,   /**< Output: Frequency (Hz) */
   REAL8 t,            /**< Time (s) at frequency */
   REAL8 m1SI,         /**< Mass of companion 1 (kg) */
   REAL8 m2SI,         /**< Mass of companion 2 (kg) */
   REAL8 chi1,         /**< Dimensionless aligned component spin 1 */
   REAL8 chi2          /**< Dimensionless aligned component spin 2 */
 );


/**
 * Compute SEOBNRv2 chirp time from interpolant assuming a single-spin.
 */
REAL8 XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(
  const REAL8 m1_SI,    /**< Mass of companion 1 [kg] */
  const REAL8 m2_SI,    /**< Mass of companion 2 [kg] */
  const REAL8 chi,      /**< Effective aligned spin */
  const REAL8 f_min     /**< Starting frequency [Hz] */
);

  
/**
 * Routine to compute the mass and spin of the final black hole given
 * the masses, spins, binding energy, and orbital angular momentum vector.
 */
int XLALSimIMRPhenSpinFinalMassSpin(REAL8 *finalMass,
				    REAL8 *finalSpin,
				    REAL8 m1,
				    REAL8 m2,
				    REAL8 s1s1,
				    REAL8 s2s2,
				    REAL8 s1L,
				    REAL8 s2L,
				    REAL8 s1s2,
				    REAL8 energy);

int XLALSimSpinInspiralGenerator(REAL8TimeSeries **hPlus,	        /**< +-polarization waveform [returned] */
				 REAL8TimeSeries **hCross,	        /**< x-polarization waveform [returned] */
				 REAL8 phi_start,                       /**< start phase */
				 REAL8 deltaT,                          /**< sampling interval */
				 REAL8 m1,                              /**< mass of companion 1 */
				 REAL8 m2,                              /**< mass of companion 2 */
				 REAL8 f_min,                           /**< start frequency */
				 REAL8 f_ref,                           /**< reference frequency */
				 REAL8 r,                               /**< distance of source */
				 REAL8 iota,                            /**< inclination of source (rad) */
				 REAL8 s1x,                             /**< x-component of dimensionless spin for object 1 */
				 REAL8 s1y,                             /**< y-component of dimensionless spin for object 1 */
				 REAL8 s1z,                             /**< z-component of dimensionless spin for object 1 */
				 REAL8 s2x,                             /**< x-component of dimensionless spin for object 2 */
				 REAL8 s2y,                             /**< y-component of dimensionless spin for object 2 */
				 REAL8 s2z,                             /**< z-component of dimensionless spin for object 2 */
				 int phaseO,                            /**< twice post-Newtonian phase order */
				 int ampO,                              /**< twice post-Newtonian amplitude order */
				 LALSimInspiralWaveformFlags *waveFlags,/**< Choice of axis for input spin params */
				 LALSimInspiralTestGRParam *testGRparams/**< Choice of axis for input spin params */
				 );

/**
 * Driver routine to compute a precessing post-Newtonian inspiral-merger-ringdown waveform
 */

int XLALSimIMRPhenSpinInspiralRDGenerator(
    REAL8TimeSeries **hplus,    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,   /**< x-polarization waveform */
    REAL8 phi0,                 /**< phase at time of peak amplitude*/
    REAL8 deltaT,               /**< sampling interval */
    REAL8 m1,                   /**< mass of companion 1 */
    REAL8 m2,                   /**< mass of companion 2 */
    REAL8 f_min,                /**< start frequency */
    REAL8 f_ref,                /**< reference frequency */
    REAL8 r,                    /**< distance of source */
    REAL8 iota,                 /**< inclination of source (rad) */
    REAL8 s1x,                  /**< x-component of dimensionless spin for object 1 */
    REAL8 s1y,                  /**< y-component of dimensionless spin for object 1 */
    REAL8 s1z,                  /**< z-component of dimensionless spin for object 1 */
    REAL8 s2x,                  /**< x-component of dimensionless spin for object 2 */
    REAL8 s2y,                  /**< y-component of dimensionless spin for object 2 */
    REAL8 s2z,                  /**< z-component of dimensionless spin for object 2 */
    int phaseO,                 /**< twice post-Newtonian phase order */
    int ampO,                   /**< twice post-Newtonian amplitude order */
    LALSimInspiralWaveformFlags *waveFlag,/**< Choice of axis for input spin params */
    LALSimInspiralTestGRParam *testGRparam  /**< Choice of axis for input spin params */
					  );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMR_H */
