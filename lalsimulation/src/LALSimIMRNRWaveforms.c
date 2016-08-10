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
UNUSED static UINT4 XLALSimInspiralNRWaveformGetSpinsFromHDF5FilePointer(
  UNUSED REAL8 *S1x,           /**< [out] Dimensionless spin1x in LAL frame */
  UNUSED REAL8 *S1y,           /**< [out] Dimensionless spin1y in LAL frame */
  UNUSED REAL8 *S1z,           /**< [out] Dimensionless spin1z in LAL frame */
  UNUSED REAL8 *S2x,           /**< [out] Dimensionless spin2x in LAL frame */
  UNUSED REAL8 *S2y,           /**< [out] Dimensionless spin2y in LAL frame */
  UNUSED REAL8 *S2z,           /**< [out] Dimensionless spin2z in LAL frame */
  UNUSED LALH5File* file       /**< Pointer to HDF5 file */
)
{
  #ifndef LAL_HDF5_ENABLED
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
  #else
  REAL8 nrSpin1x, nrSpin1y, nrSpin1z, nrSpin2x, nrSpin2y, nrSpin2z;
  REAL8 ln_hat_x, ln_hat_y, ln_hat_z, n_hat_x, n_hat_y, n_hat_z;
  XLALH5FileQueryScalarAttributeValue(&nrSpin1x, file, "spin1x");
  XLALH5FileQueryScalarAttributeValue(&nrSpin1y, file, "spin1y");
  XLALH5FileQueryScalarAttributeValue(&nrSpin1z, file, "spin1z");
  XLALH5FileQueryScalarAttributeValue(&nrSpin2x, file, "spin2x");
  XLALH5FileQueryScalarAttributeValue(&nrSpin2y, file, "spin2y");
  XLALH5FileQueryScalarAttributeValue(&nrSpin2z, file, "spin2z");

  XLALH5FileQueryScalarAttributeValue(&ln_hat_x , file, "LNhatx");
  XLALH5FileQueryScalarAttributeValue(&ln_hat_y , file, "LNhaty");
  XLALH5FileQueryScalarAttributeValue(&ln_hat_z , file, "LNhatz");

  XLALH5FileQueryScalarAttributeValue(&n_hat_x , file, "nhatx");
  XLALH5FileQueryScalarAttributeValue(&n_hat_y , file, "nhaty");
  XLALH5FileQueryScalarAttributeValue(&n_hat_z , file, "nhatz");

  *S1x = nrSpin1x * n_hat_x + nrSpin1y * n_hat_y + nrSpin1z * n_hat_z;
  *S1y = nrSpin1x * (-ln_hat_z * n_hat_y + ln_hat_y * n_hat_z)
        + nrSpin1y * (ln_hat_z * n_hat_x - ln_hat_x * n_hat_z)
        + nrSpin1z * (-ln_hat_y * n_hat_x + ln_hat_x * n_hat_y) ;
  *S1z = nrSpin1x * ln_hat_x + nrSpin1y * ln_hat_y + nrSpin1z * ln_hat_z;

  *S2x = nrSpin2x * n_hat_x + nrSpin2y * n_hat_y + nrSpin2z * n_hat_z;
  *S2y = nrSpin2x * (-ln_hat_z * n_hat_y + ln_hat_y * n_hat_z)
        + nrSpin2y * (ln_hat_z * n_hat_x - ln_hat_x * n_hat_z)
        + nrSpin2z * (-ln_hat_y * n_hat_x + ln_hat_x * n_hat_y);
  *S2z = nrSpin2x * ln_hat_x + nrSpin2y * ln_hat_y + nrSpin2z * ln_hat_z;

  return XLAL_SUCCESS;
  #endif
}

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

UNUSED static UINT4 XLALSimInspiralNRWaveformGetRotationAnglesFromH5File(
  UNUSED REAL8 *theta,            /**< Returned inclination angle of source */
  UNUSED REAL8 *psi,              /**< Returned azimuth angle of source */
  UNUSED REAL8 *calpha,           /**< Returned cosine of the polarisation angle */
  UNUSED REAL8 *salpha,           /**< Returned sine of the polarisation angle */
  UNUSED LALH5File* filepointer,  /**< Pointer to NR HDF5 file */
  UNUSED const REAL8 inclination, /**< Inclination of source */
  UNUSED const REAL8 phi_ref      /**< Orbital reference phase*/
  )
{
  #ifndef LAL_HDF5_ENABLED
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
  #else

  /* Compute the angles necessary to rotate from the intrinsic NR source frame
   * into the LAL frame. See DCC-T1600045 for details.
   */

  /* Attribute declarations go in here */
  /* Declarations */
  REAL8 orb_phase, ln_hat_x, ln_hat_y, ln_hat_z, ln_hat_norm;
  REAL8 n_hat_x, n_hat_y, n_hat_z, n_hat_norm;
  REAL8 corb_phase, sorb_phase, sinclination, cinclination;
  REAL8 ln_cross_n_x, ln_cross_n_y, ln_cross_n_z;
  REAL8 z_wave_x, z_wave_y, z_wave_z;
  REAL8 stheta, ctheta, spsi, cpsi, theta_hat_x, theta_hat_y, theta_hat_z;
  REAL8 psi_hat_x, psi_hat_y, psi_hat_z;
  REAL8 n_dot_theta, ln_cross_n_dot_theta, n_dot_psi, ln_cross_n_dot_psi;
  REAL8 y_val;

  /* Following section IV of DCC-T1600045
   * Step 1: Define Phi = phiref/2 ... I'm ignoring this, I think it is wrong
   */

  orb_phase = phi_ref;

  /* Step 2: Compute Zref
   * 2.1: Compute LN_hat from file. LN_hat = direction of orbital ang. mom.
   */

  XLALH5FileQueryScalarAttributeValue(&ln_hat_x , filepointer, "LNhatx");
  XLALH5FileQueryScalarAttributeValue(&ln_hat_y , filepointer, "LNhaty");
  XLALH5FileQueryScalarAttributeValue(&ln_hat_z , filepointer, "LNhatz");

  ln_hat_norm = sqrt(ln_hat_x * ln_hat_x + ln_hat_y * ln_hat_y + ln_hat_z * ln_hat_z);

  ln_hat_x = ln_hat_x / ln_hat_norm;
  ln_hat_y = ln_hat_y / ln_hat_norm;
  ln_hat_z = ln_hat_z / ln_hat_norm;

  /* 2.2: Compute n_hat from file.
   * n_hat = direction from object 2 to object 1
   */

  XLALH5FileQueryScalarAttributeValue(&n_hat_x , filepointer, "nhatx");
  XLALH5FileQueryScalarAttributeValue(&n_hat_y , filepointer, "nhaty");
  XLALH5FileQueryScalarAttributeValue(&n_hat_z , filepointer, "nhatz");

  n_hat_norm = sqrt(n_hat_x * n_hat_x + n_hat_y * n_hat_y + n_hat_z * n_hat_z);

  n_hat_x = n_hat_x / n_hat_norm;
  n_hat_y = n_hat_y / n_hat_norm;
  n_hat_z = n_hat_z / n_hat_norm;

  /* 2.3: Compute Z in the lal wave frame */

  corb_phase = cos(orb_phase);
  sorb_phase = sin(orb_phase);
  sinclination = sin(inclination);
  cinclination = cos(inclination);

  ln_cross_n_x = ln_hat_y * n_hat_z - ln_hat_z * n_hat_y;
  ln_cross_n_y = ln_hat_z * n_hat_x - ln_hat_x * n_hat_z;
  ln_cross_n_z = ln_hat_x * n_hat_y - ln_hat_y * n_hat_x;

  z_wave_x = sinclination * (sorb_phase * n_hat_x + corb_phase * ln_cross_n_x);
  z_wave_y = sinclination * (sorb_phase * n_hat_y + corb_phase * ln_cross_n_y);
  z_wave_z = sinclination * (sorb_phase * n_hat_z + corb_phase * ln_cross_n_z);

  z_wave_x += cinclination * ln_hat_x;
  z_wave_y += cinclination * ln_hat_y;
  z_wave_z += cinclination * ln_hat_z;

  /* Step 3.1: Extract theta and psi from Z in the lal wave frame
   * NOTE: Theta can only run between 0 and pi, so no problem with arccos here
   */

  *theta = acos(z_wave_z);

  /* Degenerate if Z_wave[2] == 1. In this case just choose psi randomly,
   * the choice will be cancelled out by alpha correction (I hope!)
   */

  if(fabs(z_wave_z - 1.0 ) < 0.000001)
  {
    *psi = 0.5;
  }
  else
  {
    /* psi can run between 0 and 2pi, but only one solution works for x and y */
    /* Possible numerical issues if z_wave_x = sin(theta) */
    if(fabs(z_wave_x / sin(*theta)) > 1.)
    {
      if(fabs(z_wave_x / sin(*theta)) < 1.00001)
      {
        if((z_wave_x * sin(*theta)) < 0.)
        {
          *psi = LAL_PI;
        }
        else
        {
          *psi = 0.;
        }
      }
    }
    else
    {
      *psi = acos(z_wave_x / sin(*theta));
    }
    y_val = sin(*psi) * sin(*theta);
    /*  If z_wave[1] is negative, flip psi so that sin(psi) goes negative
     *  while preserving cos(psi) */
    if( z_wave_y < 0.)
    {
      *psi = 2 * LAL_PI - *psi;
      y_val = sin(*psi) * sin(*theta);
    }
    if( fabs(y_val - z_wave_y) > 0.0001)
    {
      XLAL_ERROR(XLAL_EDOM, "Something's wrong in Ian's math. Tell him he's an idiot!");
    }
  }

  /* 3.2: Compute the vectors theta_hat and psi_hat */

  stheta = sin(*theta);
  ctheta = cos(*theta);
  spsi = sin(*psi);
  cpsi = cos(*psi);

  theta_hat_x = cpsi * ctheta;
  theta_hat_y = spsi * ctheta;
  theta_hat_z = - stheta;

  psi_hat_x = -spsi;
  psi_hat_y = cpsi;
  psi_hat_z = 0.0;

  /* Step 4: Compute sin(alpha) and cos(alpha) */

  n_dot_theta = n_hat_x * theta_hat_x + n_hat_y * theta_hat_y + n_hat_z * theta_hat_z;
  ln_cross_n_dot_theta = ln_cross_n_x * theta_hat_x + ln_cross_n_y * theta_hat_y + ln_cross_n_z * theta_hat_z;
  n_dot_psi = n_hat_x * psi_hat_x + n_hat_y * psi_hat_y + n_hat_z * psi_hat_z;
  ln_cross_n_dot_psi = ln_cross_n_x * psi_hat_x + ln_cross_n_y * psi_hat_y + ln_cross_n_z * psi_hat_z;

  *salpha = corb_phase * n_dot_theta - sorb_phase * ln_cross_n_dot_theta;
  *calpha = corb_phase * n_dot_psi - sorb_phase * ln_cross_n_dot_psi;

  /*  Step 5: Also useful to keep the source frame vectors as defined in
   *  equation 16 of Harald's document.
   */

  /*
   * x_source_hat[0] = corb_phase * n_hat_x - sorb_phase * ln_cross_n_x;
   * x_source_hat[1] = corb_phase * n_hat_y - sorb_phase * ln_cross_n_y;
   * x_source_hat[2] = corb_phase * n_hat_z - sorb_phase * ln_cross_n_z;
   * y_source_hat[0] = sorb_phase * n_hat_x + corb_phase * ln_cross_n_x;
   * y_source_hat[1] = sorb_phase * n_hat_y + corb_phase * ln_cross_n_y;
   * y_source_hat[2] = sorb_phase * n_hat_z + corb_phase * ln_cross_n_z;
   * z_source_hat[0] = ln_hat_x;
   * z_source_hat[1] = ln_hat_y;
   * z_source_hat[2] = ln_hat_z;
   */

  return XLAL_SUCCESS;
  #endif
}

int XLALSimInspiralNRWaveformGetSpinsFromHDF5File(
  UNUSED REAL8 *S1x,             /**< [out] Dimensionless spin1x in LAL frame */
  UNUSED REAL8 *S1y,             /**< [out] Dimensionless spin1y in LAL frame */
  UNUSED REAL8 *S1z,             /**< [out] Dimensionless spin1z in LAL frame */
  UNUSED REAL8 *S2x,             /**< [out] Dimensionless spin2x in LAL frame */
  UNUSED REAL8 *S2y,             /**< [out] Dimensionless spin2y in LAL frame */
  UNUSED REAL8 *S2z,             /**< [out] Dimensionless spin2z in LAL frame */
  UNUSED const char *NRDataFile  /**< Location of NR HDF file */
)
{
  #ifndef LAL_HDF5_ENABLED
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
  #else
  LALH5File *file;
  file = XLALH5FileOpen(NRDataFile, "r");
  XLALSimInspiralNRWaveformGetSpinsFromHDF5FilePointer(S1x, S1y, S1z, S2x, S2y,
                                                       S2z, file);

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
  REAL8 nrEta;
  REAL8 S1x, S1y, S1z, S2x, S2y, S2z;
  REAL8 Mflower, time_start_M, time_start_s, time_end_M, time_end_s;
  REAL8 chi, est_start_time, curr_h_real, curr_h_imag;
  REAL8 theta, psi, calpha, salpha;
  REAL8 distance_scale_fac;
  COMPLEX16 curr_ylm;
  REAL8TimeSeries *hplus_corr;
  REAL8TimeSeries *hcross_corr;

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

  /* Read spin metadata, L_hat, n_hat from HDF5 metadata and make sure
   * the ChooseTDWaveform() input values are consistent with the data
   * recorded in the metadata of the HDF5 file.
   * PS: This assumes that the input spins are in the LAL frame!
   */
  XLALSimInspiralNRWaveformGetSpinsFromHDF5FilePointer(&S1x, &S1y, &S1z,
                                                       &S2x, &S2y, &S2z, file);

  if (fabs(S1x - s1x) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1X IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  if (fabs(S1y - s1y) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1Y IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  if (fabs(S1z - s1z) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN1Z IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  if (fabs(S2x - s2x) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN2X IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  if (fabs(S2y - s2y) > 1E-3)
  {
     XLAL_ERROR(XLAL_EDOM, "SPIN2Y IS INCONSISTENT WITH THE NR SIMULATION.\n");
  }

  if (fabs(S2z - s2z) > 1E-3)
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

  /* Compute correct angles for hplus and hcross following LAL convention. */

  theta = psi = calpha = salpha = 0.;
  XLALSimInspiralNRWaveformGetRotationAnglesFromH5File(&theta, &psi, &calpha,
                       &salpha, file, inclination, phiRef);

  /* Create the return time series, use arbitrary epoch here. We set this
   * properly later. */
  XLALGPSAdd(&tmpEpoch, time_start_s);

  *hplus  = XLALCreateREAL8TimeSeries("H_PLUS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );

  hplus_corr = XLALCreateREAL8TimeSeries("H_PLUS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  hcross_corr = XLALCreateREAL8TimeSeries("H_CROSS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  for (curr_idx = 0; curr_idx < array_length; curr_idx++)
  {
    hplus_corr->data->data[curr_idx] = 0.0;
    hcross_corr->data->data[curr_idx] = 0.0;
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

      curr_ylm = XLALSpinWeightedSphericalHarmonic(theta, psi, -2,
                                                   model, modem);

      for (curr_idx = 0; curr_idx < array_length; curr_idx++)
      {
        curr_h_real = curr_amp->data[curr_idx]
                    * cos(curr_phase->data[curr_idx]) * distance_scale_fac;
        curr_h_imag = curr_amp->data[curr_idx]
                    * sin(curr_phase->data[curr_idx]) * distance_scale_fac;

        hplus_corr->data->data[curr_idx] = hplus_corr->data->data[curr_idx]
               + curr_h_real * creal(curr_ylm) - curr_h_imag * cimag(curr_ylm);

        hcross_corr->data->data[curr_idx] = hcross_corr->data->data[curr_idx]
               - curr_h_real * cimag(curr_ylm) - curr_h_imag * creal(curr_ylm);

      }

      XLALDestroyREAL8Vector(curr_amp);
      XLALDestroyREAL8Vector(curr_phase);

    }

  }

 /* Correct for the "alpha" angle as given in T1600045 to translate
  * from the NR wave frame to LAL wave-frame
  * Helper time series needed.
  */

  for (curr_idx = 0; curr_idx < array_length; curr_idx++)
  {
    (*hplus)->data->data[curr_idx] =
          (calpha*calpha - salpha*salpha) * hplus_corr->data->data[curr_idx]
          - 2.0*calpha*salpha * hcross_corr->data->data[curr_idx];

    (*hcross)->data->data[curr_idx] =
          + 2.0*calpha*salpha * hplus_corr->data->data[curr_idx]
        + (calpha*calpha - salpha*salpha) * hcross_corr->data->data[curr_idx];
  }

  XLALDestroyREAL8TimeSeries(hplus_corr);
  XLALDestroyREAL8TimeSeries(hcross_corr);
  XLALH5FileClose(file);

  return XLAL_SUCCESS;
  #endif
}
