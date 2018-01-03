#ifndef _LALSIM_IMR_LACKEY_TIDAL_2013_H
#define _LALSIM_IMR_LACKEY_TIDAL_2013_H

static void tidalPNAmplitudeCoefficient(
  double *C,
  const double eta,
  const double chi_BH,
  const double Lambda
);

static double tidalCorrectionAmplitude(
  const double mf,
  const double C,
  const double eta,
  const double Lambda
);

// precompute a0, a1 and G which do not depend on frequency
static void tidalPNPhaseCoefficients(
  double *a0,
  double *a1,
  double *G,
  const double eta,
  const double chi_BH,
  const double Lambda
);

static double tidalPNPhase(
  const double mf,
  const double a0,
  const double a1,
  const double eta
);

static double tidalPNPhaseDeriv(
  const double mf,
  const double a0,
  const double a1,
  const double eta
);

// Implements Eq. 34 of Lackey et al
static double tidalCorrectionPhase(
  const double mf,
  const double a0,
  const double a1,
  const double G,
  const double eta,
  const double Lambda
);

int LackeyTidal2013SEOBNRv2ROMCore(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 mBH_SI,                                 /**< Mass of black hole (kg) */
  REAL8 mNS_SI,                                 /**< Mass of neutron star (kg) */
  REAL8 chi_BH,                                 /**< Dimensionless aligned component spin of the BH */
  REAL8 Lambda,                                 /**< Dimensionless tidal deformability (Eq 1  of Lackey et al) */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF                                  /**< Sampling frequency (Hz) */
);

#endif /* _LALSIM_IMR_LACKEY_TIDAL_2013_H */