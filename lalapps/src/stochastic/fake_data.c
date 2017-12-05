/*
 * fake_data.c - SGWB Standalone Analysis Pipeline
 *             - Fake Data Functions
 *
 * Copyright (C) 2002-2006 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#include "fake_data.h"
#include "sgwb.h"


/* externally declared variables */
extern int vrbflg;

/* generate random noise */
REAL4TimeSeries *generate_random_noise(LALStatus *status,
    INT4 duration,
    INT4 sample_rate)
{
  /* variables */
  REAL4TimeSeries *series;
  LIGOTimeGPS start;
  REAL8 delta_t;
  RandomParams *noise_params = NULL;
  struct timeval tv;
#if 0
  UINT4 i;
#endif

  /* initialise time */
  start.gpsSeconds = 0;
  start.gpsNanoSeconds = 0;

  /* get deltaT from sample_rate */
  delta_t = 1.0/(REAL8)(sample_rate);

  if (vrbflg)
    fprintf(stdout, "Allocating memory for random noise...\n");

  /* create and initialise time series */
  series = XLALCreateREAL4TimeSeries("noise", &start, 0, delta_t, \
      &lalADCCountUnit, duration * sample_rate);

  /* get current time, for random seed */
  gettimeofday(&tv, NULL);

  if (vrbflg)
    fprintf(stdout, "Generating random noise...\n");

  /* generate noise_params */
  LAL_CALL(LALCreateRandomParams(status, &noise_params, tv.tv_usec), status);

  /* generate random noise */
  LAL_CALL(LALNormalDeviates(status, series->data, noise_params), status);

  /* destroy noise_params */
  LAL_CALL(LALDestroyRandomParams(status, &noise_params), status);

  /* scale random noise */
#if 0
  for (i = 0; i < series->data->length; i++)
    series->data->data[i] *= 1e-18;
#endif

  return(series);
}

/* generate fake detector output */
SSSimStochBGOutput *generate_fake_detector_output(LALStatus *status,
    REAL4TimeSeries *noise_one,
    REAL4TimeSeries *noise_two,
    REAL8 deltaF,
    REAL8 f_min,
    REAL8 f_max)
{
  /* variables */
  REAL4TimeSeries *series_one;
  REAL4TimeSeries *series_two;
  SSSimStochBGInput input;
  SSSimStochBGParams params;
  SSSimStochBGOutput *output = NULL;
  COMPLEX8FrequencySeries *response_one = NULL;
  COMPLEX8FrequencySeries *response_two = NULL;
  LIGOTimeGPS start;
  INT4 freq_length;
  struct timeval tv;
  UINT4 i;
  REAL4FrequencySeries *omega;

  /* initialise epoch */
  start.gpsSeconds = 0;
  start.gpsNanoSeconds = 0;

  /* get frequency length */
  freq_length = (noise_one->data->length / 2) + 1;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for fake detector output...\n");

  /* create and initialise time series' */
  series_one = XLALCreateREAL4TimeSeries("one", &noise_one->epoch, \
      noise_one->f0, noise_one->deltaT, &noise_one->sampleUnits, \
      noise_one->data->length);
  series_two = XLALCreateREAL4TimeSeries("two", &noise_two->epoch, \
      noise_two->f0, noise_two->deltaT, &noise_two->sampleUnits, \
      noise_two->data->length);

  if (vrbflg)
    fprintf(stdout, "Generating unity response functions...\n");

  /* generate unity response functions */
  response_one = unity_response(start, 0, deltaF, lalDimensionlessUnit, \
      freq_length);
  response_two = unity_response(start, 0, deltaF, lalDimensionlessUnit, \
      freq_length);

  /* generate omega */
  if (vrbflg)
    fprintf(stdout, "Generating spectrum for optimal filter...\n");
  omega = omega_gw(status, 0, (f_max - f_min) / 2, 1, freq_length, \
      0, deltaF);

  /* get current time, for random seed */
  gettimeofday(&tv, NULL);

  /* setup params for fake signal generation */
  params.length = noise_one->data->length;
  params.deltaT = noise_one->deltaT;
  params.seed = tv.tv_usec;
  params.detectorOne = lalCachedDetectors[LALDetectorIndexLHODIFF];
  params.detectorTwo = lalCachedDetectors[LALDetectorIndexLHODIFF];
  params.SSimStochBGTimeSeries1Unit = lalStrainUnit;
  params.SSimStochBGTimeSeries2Unit = lalStrainUnit;

  /* setup intputs for fake signal generation */
  input.omegaGW = omega;
  input.whiteningFilter1 = response_one;
  input.whiteningFilter2 = response_two;

  /* setup output for fake signal generation */
  output->SSimStochBG1 = series_one;
  output->SSimStochBG2 = series_two;

  if (vrbflg)
    fprintf(stdout, "Generating fake signal...\n");

  /* generate fake signal */
  LAL_CALL(LALSSSimStochBGTimeSeries(status, output, &input, &params), \
      status);

  if (vrbflg)
    fprintf(stdout, "Injecting fake signal into random noise...\n");

  /* inject signal into noise */
  for (i = 0; i < series_one->data->length; i++)
  {
    series_one->data->data[i] += noise_one->data->data[i];
    series_two->data->data[i] += noise_two->data->data[i];
  }

  return(output);
}

/*
 * vim: et
 */
