/*
 * Copyright (C) 2012 Prayush Kumar, Frank Ohme
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


/* The paper refered to here as the Main paper, is Phys. Rev. D 82, 064016 (2010)
 * */

#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/StringInput.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>

#include "LALSimIMRPhenomC_internals.c"


/**
 * private function prototypes; all internal functions use solar masses.
 *
 */

static int IMRPhenomCGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 deltaF, const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomCParams *params);

/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomC in the frequency domain.
 *
 * Reference: http://arxiv.org/pdf/1005.3306v3.pdf
 * - Waveform: Eq.(5.3)-(5.13)
 * - Coefficients: Eq.(5.14) and Table II
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomCGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< orbital phase at peak (rad) */
    const REAL8 deltaF,                /**< sampling interval (Hz) */
    const REAL8 m1_SI,                 /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< mass of companion 2 (kg) */
    const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< end frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 distance               /**< distance of source (m) */
) {
  BBHPhenomCParams *params;
  int status;
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* If spins are above 0.9 or below -0.9, throw an error */
  if (chi > 0.9 || chi < -0.9){
      XLALPrintError("Spins outside the range [-0.9,0.9] are not supported\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* If mass ratio is above 4 and below 20, give a warning, and if it is above
   * 20, throw an error */
  REAL8 q = m1/m2;
  if (q > 1.0)
    q = 1./q;

  if (q > 20.0){
      XLALPrintError("Mass ratio is way outside the calibration range. m1/m2 should be <= 20.\n");
      XLAL_ERROR(XLAL_EDOM);
  }
  else if (q > 4.0){
      XLALPrintWarning("Warning: The model is only calibrated for m1/m2 <= 4.\n");
  }

  /* phenomenological parameters*/
  params = ComputeIMRPhenomCParams(m1, m2, chi);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("(fCut = 0.15M) <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  f_max_prime = (f_max_prime > params->fCut) ? params->fCut : f_max_prime;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  status = IMRPhenomCGenerateFD(htilde, phi0, deltaF, m1, m2, f_min, f_max_prime, distance, params);
  LALFree(params);
  return status;
}

/*
 * Convenience function to quickly find the default final
 * frequency.
 */
double XLALSimIMRPhenomCGetFinalFreq(
    const REAL8 m1,
    const REAL8 m2,
    const REAL8 chi
) {
    BBHPhenomCParams *phenomParams;
    phenomParams = ComputeIMRPhenomCParams(m1, m2, chi);
    return phenomParams->fCut;
}

/***********************************************************************************/
/* The following private function generates IMRPhenomC frequency-domain waveforms  */
/* given coefficients */
/***********************************************************************************/

static int IMRPhenomCGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at peak */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    //const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance to source (m) */
    const BBHPhenomCParams *params      /**< from ComputeIMRPhenomCParams */
) {
  static LIGOTimeGPS ligotimegps_zero = {0, 0};
  size_t i;
  INT4 errcode;

  const REAL8 M = m1 + m2;
  const REAL8 eta = m1 * m2 / (M * M);

  /* Memory to temporarily store components of amplitude and phase */
  /*
  REAL8 phSPA, phPM, phRD, aPM, aRD;
  REAL8 wPlusf1, wPlusf2, wMinusf1, wMinusf2, wPlusf0, wMinusf0;
  */
  REAL8 phPhenomC = 0.0;
  REAL8 aPhenomC = 0.0;

  /* compute the amplitude pre-factor */
  REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

  /* allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0,
      deltaF, &lalStrainUnit, n);
  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

  /* now generate the waveform */
  size_t ind_max = (size_t) (f_max / deltaF);
  for (i = (size_t) (f_min / deltaF); i < ind_max; i++)
  {
    REAL8 f = i * deltaF;

    errcode = IMRPhenomCGenerateAmpPhase( &aPhenomC, &phPhenomC, f, eta, params );
    if( errcode != XLAL_SUCCESS )
      XLAL_ERROR(XLAL_EFUNC);
    phPhenomC -= 2.*phi0; // factor of 2 b/c phi0 is orbital phase

    /* generate the waveform */
    ((*htilde)->data->data)[i] = amp0 * aPhenomC * cos(phPhenomC);
    ((*htilde)->data->data)[i] += -I * amp0 * aPhenomC * sin(phPhenomC);
  }

  return XLAL_SUCCESS;
}
