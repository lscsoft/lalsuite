#ifndef _LALSIM_IMR_SEOBNRv4_ROM_NRTidal_H
#define _LALSIM_IMR_SEOBNRv4_ROM_NRTidal_H

static int SEOBNRv4ROM_NRTidal_Core(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1_SI,                                  /**< Mass of neutron star 1 (kg) */
  REAL8 m2_SI,                                  /**< Mass of neutron star 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin of NS 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin of NS 2 */
  REAL8 Lambda1,                                /**< Dimensionless tidal deformability of NS 1 */
  REAL8 Lambda2,                                /**< Dimensionless tidal deformability of NS 2 */
  const REAL8Sequence *freqs_in,                /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF                                  /**< Sampling frequency (Hz) */
);

#endif
