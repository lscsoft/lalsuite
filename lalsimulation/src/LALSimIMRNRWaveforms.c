/*
* Copyright (C) 2015 Ian Harry, Patricia Schmidt
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <libgen.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConfig.h>
#include <lal/SphericalHarmonics.h>

#include <lal/H5FileIO.h>

#include "LALSimIMRSEOBNRROMUtilities.c"

/* Everything needs to be declared as unused in case HDF is not enabled. */
UNUSED static UINT4 XLALSimInspiralNRWaveformGetDataFromHDF5File(
  UNUSED REAL8Vector** output,            /**< Returned vector uncompressed */
  UNUSED LALH5File* pointer,              /**< Pointer to HDF5 file */
  UNUSED REAL8 totalMass,                 /**< Total mass of system for scaling */
  UNUSED REAL8 startTime,                 /**< Start time of veturn vector */
  UNUSED size_t length,                   /**< Length of returned vector */
  UNUSED REAL8 deltaT,                    /**< Sample rate of returned vector */
  UNUSED const char *keyName              /**< Name of vector to uncompress */
  )
{
  #ifndef LAL_HDF5_ENABLED
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
  #else
  UINT4 idx;
  size_t comp_data_length;
  REAL8 massTime;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  gsl_vector *knotsVector, *dataVector;
  LALH5File *group = XLALH5GroupOpen(pointer, keyName);
  knotsVector=dataVector=NULL;

  ReadHDF5RealVectorDataset(group, "knots", &knotsVector);
  ReadHDF5RealVectorDataset(group, "data", &dataVector);

  *output = XLALCreateREAL8Vector(length);

  comp_data_length = dataVector->size;
  /* SPLINE STUFF */
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, comp_data_length);
  gsl_spline_init(spline, knotsVector->data, dataVector->data,
                  comp_data_length);

  for (idx = 0; idx < length; idx++)
  {
    massTime = (startTime + idx*deltaT) / (totalMass * LAL_MTSUN_SI);
    (*output)->data[idx] = gsl_spline_eval(spline, massTime, acc);
  }

  gsl_vector_free(knotsVector);
  gsl_vector_free(dataVector);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return XLAL_SUCCESS;
  #endif
}

/* Everything needs to be declared as unused in case HDF is not enabled. */
int XLALSimInspiralNRWaveformGetHplusHcross(
        UNUSED REAL8TimeSeries **hplus,        /**< Output h_+ vector */
        UNUSED REAL8TimeSeries **hcross,       /**< Output h_x vector */
        UNUSED REAL8 phiRef,                   /**< orbital phase at reference pt. */
        UNUSED REAL8 inclination,              /**< inclination angle */
        UNUSED REAL8 deltaT,                   /**< sampling interval (s) */
        UNUSED REAL8 m1,                       /**< mass of companion 1 (kg) */
        UNUSED REAL8 m2,                       /**< mass of companion 2 (kg) */
        UNUSED REAL8 r,                        /**< distance of source (m) */
        UNUSED REAL8 fStart,                   /**< start GW frequency (Hz) */
        UNUSED REAL8 fRef,                     /**< reference GW frequency (Hz) */
        UNUSED REAL8 s1x,                      /**< initial value of S1x */
        UNUSED REAL8 s1y,                      /**< initial value of S1y */
        UNUSED REAL8 s1z,                      /**< initial value of S1z */
        UNUSED REAL8 s2x,                      /**< initial value of S2x */
        UNUSED REAL8 s2y,                      /**< initial value of S2y */
        UNUSED REAL8 s2z,                      /**< initial value of S2z */
        UNUSED const char *NRDataFile          /**< Location of NR HDF file */
        )
{
  #ifndef LAL_HDF5_ENABLED
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
  #else
  /* Declarations */
  UINT4 curr_idx;
  INT4 model, modem;
  size_t array_length;
  REAL8 nrEta, nrSpin1x, nrSpin1y, nrSpin1z, nrSpin2x, nrSpin2y, nrSpin2z;
  REAL8 Mflower, time_start_M, time_start_s, time_end_M, time_end_s;
  REAL8 chi, est_start_time, nrPhiRef, phi, curr_h_real, curr_h_imag;
  REAL8 distance_scale_fac;
  COMPLEX16 curr_ylm;
  /* These keys follow a strict formulation and cannot be longer than 11
   * characters */
  char amp_key[20];
  char phase_key[20];
  gsl_vector *tmpVector=NULL;
  LALH5File *file, *group;
  LIGOTimeGPS tmpEpoch = LIGOTIMEGPSZERO;
  REAL8Vector *curr_amp, *curr_phase;

  /* Use solar masses for units. NR files will use
   * solar masses as well, so easier for that conversion
   */
  m1 = m1 / LAL_MSUN_SI;
  m2 = m2 / LAL_MSUN_SI;

  file = XLALH5FileOpen(NRDataFile, "r");

  /* Sanity checks on physical parameters passed to waveform
   * generator to guarantee consistency with NR data file.
   */
  XLALH5FileQueryScalarAttributeValue(&nrEta, file, "eta");
  if (fabs((m1 * m2) / pow((m1 + m2),2.0) - nrEta) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "MASSES ARE INCONSISTENT WITH THE MASS RATIO OF THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin1x, file, "spin1x");
  if (fabs(nrSpin1x - s1x) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1X IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin1y, file, "spin1y");
  if (fabs(nrSpin1y - s1y) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1Y IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin1z, file, "spin1z");
  if (fabs(nrSpin1z - s1z) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1Z IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin2x, file, "spin2x");
  if (fabs(nrSpin2x - s2x) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN2X IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin2y, file, "spin2y");
  if (fabs(nrSpin2y - s2y) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN2Y IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  XLALH5FileQueryScalarAttributeValue(&nrSpin2z, file, "spin2z");
  if (fabs(nrSpin2z - s2z) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN2Z IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  /* First estimate the length of time series that is needed.
   * Demand that 22 mode that is present and use that to figure this out
   */
  XLALH5FileQueryScalarAttributeValue(&Mflower, file, "f_lower_at_1MSUN");
  /* Figure out start time of data */
  group = XLALH5GroupOpen(file, "amp_l2_m2");
  ReadHDF5RealVectorDataset(group, "knots", &tmpVector);
  time_start_M = (REAL8)(gsl_vector_get(tmpVector, 0));
  time_end_M = (REAL8)(gsl_vector_get(tmpVector, tmpVector->size - 1));
  gsl_vector_free(tmpVector);
  time_start_s = time_start_M * (m1 + m2) * LAL_MTSUN_SI;
  time_end_s = time_end_M * (m1 + m2) * LAL_MTSUN_SI;

  /* We don't want to return the *entire* waveform if it will be much longer
   * than the specified f_lower. Therefore guess waveform length using
   * the SEOBNR_ROM function and add 10% for safety.
   * FIXME: Is this correct for precessing waveforms?
   */
  chi = XLALSimIMRPhenomBComputeChi(m1, m2, s1z, s2z);
  est_start_time = XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(m1 * LAL_MSUN_SI,
                                                m2 * LAL_MSUN_SI, chi, fStart);
  est_start_time = (-est_start_time) * 1.1;
  if (est_start_time > time_start_s)
  {
    /* Restrict start time of waveform */
    time_start_s = est_start_time;
    time_start_M = time_start_s / ((m1 + m2) * LAL_MTSUN_SI);
  }
  else if (fStart < Mflower / (m1 + m2) )
  {
     XLAL_ERROR(XLAL_EDOM, "WAVEFORM IS NOT LONG ENOUGH TO REACH f_low. %e %e %e",
                fStart, Mflower, Mflower / (m1 + m2));
  }

  array_length = (UINT4)(ceil( (time_end_s - time_start_s) / deltaT));

  /* Create the return time series, use arbitrary epoch here. We set this
   * properly later. */
  XLALGPSAdd(&tmpEpoch, time_start_s);
  *hplus  = XLALCreateREAL8TimeSeries("H_PLUS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  for (curr_idx = 0; curr_idx < array_length; curr_idx++)
  {
    (*hplus)->data->data[curr_idx] = 0.;
    (*hcross)->data->data[curr_idx] = 0.;
  }

  /* Create the distance scale factor */
  distance_scale_fac = (m1 + m2) * LAL_MRSUN_SI / r;

  /* Generate the waveform */
  for (model=2; model < 9 ; model++)
  {
    for (modem=-model; modem < (model+1); modem++)
    {
      snprintf(amp_key, sizeof(amp_key), "amp_l%d_m%d", model, modem);
      snprintf(phase_key, sizeof(phase_key), "phase_l%d_m%d", model, modem);

      /* Check that both groups exist */
      if (XLALH5CheckGroupExists(file, amp_key) == 0)
      {
        continue;
      }
      if (XLALH5CheckGroupExists(file, phase_key) == 0)
      {
        continue;
      }

      /* Get amplitude and phase from file */
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_amp, file, (m1 + m2),
                                  time_start_s, array_length, deltaT, amp_key);
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_phase, file, (m1 + m2),
                                time_start_s, array_length, deltaT, phase_key);

      XLALH5FileQueryScalarAttributeValue(&nrPhiRef, file, "coa_phase");
      phi = nrPhiRef + phiRef;
      curr_ylm = XLALSpinWeightedSphericalHarmonic(inclination, phi, -2,
                                                   model, modem);
      for (curr_idx = 0; curr_idx < array_length; curr_idx++)
      {
        curr_h_real = curr_amp->data[curr_idx]
                    * cos(curr_phase->data[curr_idx]) * distance_scale_fac;
        curr_h_imag = curr_amp->data[curr_idx]
                    * sin(curr_phase->data[curr_idx]) * distance_scale_fac;
        (*hplus)->data->data[curr_idx] = (*hplus)->data->data[curr_idx]
               + curr_h_real * creal(curr_ylm) + curr_h_imag * cimag(curr_ylm);
        (*hcross)->data->data[curr_idx] = (*hcross)->data->data[curr_idx]
               + curr_h_real * cimag(curr_ylm) - curr_h_imag * creal(curr_ylm);
      }
      XLALDestroyREAL8Vector(curr_amp);
      XLALDestroyREAL8Vector(curr_phase);
    }
  }

  XLALH5FileClose(file);

  return XLAL_SUCCESS;
  #endif
}
