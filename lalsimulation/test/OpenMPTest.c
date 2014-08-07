/*
 *  Copyright (C) 2014 Leo Singer, Michael Puerrer
 *
 *  Check that OMP enable waveforms give identical answers no matter how
 *  many OpenMP threads are used.
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


#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Units.h>

#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* Return 1 if two COMPLEX16FrequencySeries differ,
 * or 0 if they are identical. Print a message to stderr
 * describing the first field that differs. */
static int series_differ(
    COMPLEX16FrequencySeries *a,
    COMPLEX16FrequencySeries *b
) {
    int ret = 1;

    if (a==NULL && b==NULL) // Catch unallocated frequency series
      return 0;

    if (strcmp(a->name, b->name) != 0)
        fputs("name differs", stderr);
    else if (XLALGPSCmp(&a->epoch, &b->epoch) != 0)
        fputs("epoch differs", stderr);
    else if (a->f0 != b->f0)
        fputs("f0 differs", stderr);
    else if (a->deltaF != b->deltaF)
        fputs("deltaF differs", stderr);
    else if (XLALUnitCompare(&a->sampleUnits, &b->sampleUnits) != 0)
        fputs("sampleUnits differs", stderr);
    else if (a->data->length != b->data->length)
        fputs("length differs", stderr);
    else if (memcmp(a->data->data, b->data->data,
                    a->data->length * sizeof(a->data->data[0])) != 0)
        fputs("data differs", stderr);
    else
        ret = 0;

    return ret;
}

typedef enum {
  OMP_FIRST = 0,
  OMP_TaylorF2,
  OMP_IMRPhenomP,
  OMP_LAST
} OpenMPCapableWaveforms;

static void GenerateOMPWaveform(OpenMPCapableWaveforms wf, COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde) {
  REAL8 m1 = 1.4;
  REAL8 m2 = 5.6;
  REAL8 m1_SI = m1 * LAL_MSUN_SI;
  REAL8 m2_SI = m2 * LAL_MSUN_SI;
  REAL8 s1x = 0.3;
  REAL8 s1y = 0;
  REAL8 s1z = 0.45;
  REAL8 s2x = 0;
  REAL8 s2y = 0;
  REAL8 s2z = 0.45;
  REAL8 lnhatx = sin(0.4);
  REAL8 lnhaty = 0;
  REAL8 lnhatz = cos(0.4);
  REAL8 f_min = 10;
  REAL8 f_ref = 10;
  REAL8 f_max = 0;
  REAL8 phi_ref = 0;
  REAL8 deltaF = 1. / 1024;
  REAL8 distance = 1e6 * LAL_PC_SI;
  REAL8 lambda1 = 0;
  REAL8 lambda2 = 0;
  LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_0PN;
  LALSimInspiralTidalOrder tidalO = LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN;
  INT4 phaseO = LAL_PNORDER_THREE_POINT_FIVE;
  INT4 amplitudeO = LAL_PNORDER_THREE_POINT_FIVE;

  REAL8 eta, chi_eff, chip, thetaJ, phiJ, alpha0;

  switch (wf) {
    case OMP_TaylorF2:
      XLALSimInspiralTaylorF2(
          &hptilde, phi_ref, deltaF,
          m1_SI, m2_SI,
          s1z, s2z, f_min, f_max, f_ref, distance,
          lambda1, lambda2,
          spinO,
          tidalO,
          phaseO,
          amplitudeO);
      break;
    case OMP_IMRPhenomP:
      XLALSimIMRPhenomPCalculateModelParameters(
          &chi_eff,           /**< Output: Effective aligned spin */
          &chip,              /**< Output: Effective spin in the orbital plane */
          &eta,               /**< Output: Symmetric mass-ratio */
          &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
          &phiJ,              /**< Output: Angle of J0 in the plane of the sky */
          &alpha0,            /**< Output: Initial value of alpha angle */
          m1_SI,              /**< Mass of companion 1 (kg) */
          m2_SI,              /**< Mass of companion 2 (kg) */
          f_ref,              /**< Starting GW frequency (Hz) */
          lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
          lnhaty,             /**< Initial value of LNhaty */
          lnhatz,             /**< Initial value of LNhatz */
          s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
          s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
          s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
          s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
          s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
          s2z);               /**< Initial value of s2z: dimensionless spin of larger BH */
        XLALSimIMRPhenomP(
          &hptilde,           /**< Frequency-domain waveform h+ */
          &hctilde,           /**< Frequency-domain waveform hx */
          chi_eff,            /**< Effective aligned spin */
          chip,               /**< Effective spin in the orbital plane */
          eta,                /**< Symmetric mass-ratio */
          thetaJ,             /**< Angle between J0 and line of sight (z-direction) */
          phiJ,               /**< Angle of J0 in the plane of the sky */
          m1_SI + m2_SI,      /**< Total mass of binary (kg) */
          distance,           /**< Distance of source (m) */
          alpha0,             /**< Initial value of alpha angle */
          phi_ref,            /**< Orbital coalescence phase (rad) */
          deltaF,             /**< Sampling frequency (Hz) */
          f_min,              /**< Starting GW frequency (Hz) */
          f_max,              /**< End frequency; 0 defaults to ringdown cutoff freq */
          f_ref);             /**< Reference frequency */
      break;
    default:
      XLALPrintError("Error: waveform %d not listed under OpenMPCapableWaveforms.\n", wf);
  }

}

int main (int argc, char **argv) {
    int num_threads;
    COMPLEX16FrequencySeries *base_hptilde = NULL;
    COMPLEX16FrequencySeries *base_hctilde = NULL;


    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;


    /* Loop over all OMP capable waveforms we know */
    for (OpenMPCapableWaveforms wf = OMP_FIRST+1; wf < OMP_LAST; wf++)
    {

      /* Check that using 2-8 threads gives an answer that is identical to using 1 thread. */
      for (num_threads = 1; num_threads < 8; num_threads++)
      {
	  COMPLEX16FrequencySeries *hptilde = NULL;
	  COMPLEX16FrequencySeries *hctilde = NULL;

	  omp_set_num_threads(num_threads);

	  GenerateOMPWaveform(wf, hptilde, hctilde);

	  if (num_threads == 1) {
	      base_hptilde = hptilde;
	      base_hctilde = hctilde;
	  }
	  else if (series_differ(base_hptilde, hptilde) || series_differ(base_hctilde, hctilde)) {
	      XLALPrintError("Error: frequency series differ for waveform %d.\n", wf);
	      return 1;
	  }
	  else {
	      XLALDestroyCOMPLEX16FrequencySeries(hptilde);
	      XLALDestroyCOMPLEX16FrequencySeries(hctilde);
	  }
      }
   }

   return 0;
}
