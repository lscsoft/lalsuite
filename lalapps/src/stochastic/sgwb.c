/*
 * sgwb.c - SGWB Standalone Analysis Pipeline
 *        - Main Stochastic Search Functions
 * 
 * Copyright (C) 2002-2006 Adam Mercer
 * Copyright (C) 2003-2004 Tania Regimbau
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
 * $Id$
 */

#include "sgwb.h"
#include "misc.h"
#include "data_output.h"

NRCSID(DATAINPUTC, "$Id$");
RCSID("$Id$");

/* externally declared variables */
extern int overlap_hann_flag;
extern int debug_flag;
extern int cc_spectra_flag;
extern int middle_segment_flag;
extern int vrbflg;

/* generate a data window */
REAL4Window *data_window(REAL8 delta_t,
    INT4 length,
    INT4 hann_duration)
{
  /* variables */
  REAL4Window *window = NULL;
  REAL4Window *hann = NULL;
  INT4 hann_length;
  INT4 i;

  /* get length of hann segment requested */
  hann_length = (INT4)(hann_duration / delta_t);

  if (hann_length == 0)
  {
    /* rectangular window requested */
    window = XLALCreateRectangularREAL4Window(length);
  }
  else if (hann_length == length)
  {
    /* pure hann window requested */
    window = XLALCreateHannREAL4Window(length);
  }
  else if ((hann_length > 0) && (hann_length < length))
  {
    window = XLALCreateRectangularREAL4Window(length);
    hann =  XLALCreateHannREAL4Window(hann_length);

    /* construct tukey window */
    for (i = 0; i < hann_length / 2; i++)
      window->data->data[i] = hann->data->data[i];
    for (i = hann_length / 2; i < hann_length; i++)
      window->data->data[length - hann_length + i] = hann->data->data[i];

    /* free memory for hann window */
    XLALDestroyREAL4Window(hann);
  }
  else
  {
    fprintf(stderr, "Invalid hann_length to data_window()...\n");
    exit(1);
  }

  return(window);
}

/* return the spectrum */
REAL4FrequencySeries *omega_gw(LALStatus *status,
    REAL4 exponent,
    REAL8 f_ref,
    REAL4 omega_ref,
    UINT4 length,
    REAL8 f0,
    REAL8 delta_f)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOmegaGWParameters omega_params;
  LIGOTimeGPS epoch;

  /* set epoch */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* create and initialise frequency series */
  series = XLALCreateREAL4FrequencySeries("OmegaGW", &epoch, f0, delta_f, \
      &lalDimensionlessUnit, length);

  /* set parameters */
  omega_params.alpha = exponent;
  omega_params.fRef = f_ref;
  omega_params.omegaRef = omega_ref;
  omega_params.length = length;
  omega_params.f0 = f0;
  omega_params.deltaF = delta_f;

  /* calculate spectrum */
  LAL_CALL(LALStochasticOmegaGW(status, series, &omega_params), status);

  return(series);
}

/* return the overlap reduction function */
REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length,
    REAL8 f0,
    REAL8 delta_f,
    INT4 site_one,
    INT4 site_two)
{
  /* variables */
  REAL4FrequencySeries *series;
  OverlapReductionFunctionParameters overlap_params;
  LALDetectorPair detectors;
  LIGOTimeGPS epoch;

  /* set epoch */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* create and initialise frequency series */
  series = XLALCreateREAL4FrequencySeries("Overlap", &epoch, f0, delta_f, \
      &lalDimensionlessUnit, length);

  /* set parameters */
  overlap_params.length = length;
  overlap_params.f0 = f0;
  overlap_params.deltaF = delta_f;

  /* set detectors */
  detectors.detectorOne = lalCachedDetectors[site_one];
  detectors.detectorTwo = lalCachedDetectors[site_two];

  /* calculate overlap reduction function */
  LAL_CALL(LALOverlapReductionFunction(status, series, &detectors, \
        &overlap_params), status);

  return(series);
}

/* calculate and return the inverse noise */
REAL4FrequencySeries *inverse_noise(LALStatus *status,
    REAL4FrequencySeries *psd,
    COMPLEX8FrequencySeries *response)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticInverseNoiseInput input;
  StochasticInverseNoiseCalOutput output;

  /* allocate memory */
  series = XLALCreateREAL4FrequencySeries("calPSD", &response->epoch, \
      response->f0, response->deltaF, &lalDimensionlessUnit, \
      response->data->length);

  /* set input */
  input.unCalibratedNoisePSD = psd;
  input.responseFunction = response;

  /* set output */
  output.calibratedInverseNoisePSD = series;

  /* generate inverse noise */
  LAL_CALL(LALStochasticInverseNoiseCal(status, &output, &input), status);

  return(series);
}

/* calculate and return the optimal filter */
REAL4FrequencySeries *optimal_filter(LALStatus *status,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4FrequencySeries *psd_one,
    REAL4FrequencySeries *psd_two,
    REAL4Window *window,
    REAL8 *sigma,
    REAL8 f_ref,
    INT4 segment_duration)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOptimalFilterNormalizationInput norm_input;
  StochasticOptimalFilterNormalizationOutput norm_output;
  StochasticOptimalFilterNormalizationParameters norm_params;
  StochasticOptimalFilterCalInput input;
  REAL4WithUnits normalisation;
  REAL4WithUnits variance;

  /* set parameters for normalisation */
  norm_params.fRef = f_ref;
  norm_params.heterodyned = 0;
  norm_params.window1 = window->data;
  norm_params.window2 = window->data;

  /* set inputs for normalisation */
  norm_input.overlapReductionFunction = overlap;
  norm_input.omegaGW = omega;
  norm_input.inverseNoisePSD1 = psd_one;
  norm_input.inverseNoisePSD2 = psd_two;

  /* set normalisation output */
  norm_output.normalization = &normalisation;
  norm_output.variance = &variance;

  /* calculate variance and normalisation for the optimal filter */
  LAL_CALL(LALStochasticOptimalFilterNormalization(status, \
        &norm_output, &norm_input, &norm_params), status);

  /* get theoretical sigma */
  *sigma = sqrt((REAL8)(segment_duration * variance.value * \
        pow(10, variance.units.powerOfTen)));

  /* allocate memory */
  series = XLALCreateREAL4FrequencySeries("filter", &psd_one->epoch, \
      psd_one->f0, psd_one->deltaF, &lalDimensionlessUnit, \
      psd_one->data->length);

  /* set input */
  input.overlapReductionFunction = overlap;
  input.omegaGW = omega;
  input.calibratedInverseNoisePSD1 = psd_one;
  input.calibratedInverseNoisePSD2 = psd_two;

  /* generate optimal filter */
  LAL_CALL(LALStochasticOptimalFilterCal(status, series, &input, \
        &normalisation), status);

  return(series);
}

/* estimate the PSD */
REAL4FrequencySeries *estimate_psd(REAL4TimeSeries *series,
    REAL8 f0,
    INT4 shrink_length,
    INT4 psd_window_duration)
{
  /* variables */
  REAL4FrequencySeries *psd;
  REAL4Window *window;
  RealFFTPlan *plan = NULL;
  UINT4 length;
  UINT4 overlap;
  REAL8 delta_f;
  UINT4 psd_length;

  /* lengths */
  length = psd_window_duration / series->deltaT;
  overlap = length / 2;
  delta_f = 1./(REAL8)psd_window_duration;
  psd_length = (length / 2) + 1;

  /* allocate memory */
  psd = XLALCreateREAL4FrequencySeries("psd", &series->epoch, \
      series->f0, delta_f, &lalDimensionlessUnit, psd_length);

  /* create window for PSD estimation */
  window = XLALCreateHannREAL4Window(length);

  /* create fft plan for PSD estimation */
  plan = XLALCreateForwardREAL4FFTPlan(length, 0);

  /* esimate PSD */
  XLALREAL4AverageSpectrumWelch(psd, series, length, overlap, window, plan);

  /* destroy fft plan */
  XLALDestroyREAL4FFTPlan(plan);

  /* free memory for window */
  XLALDestroyREAL4Window(window);

  /* reduce to relevant frequency range */
  XLALShrinkREAL4FrequencySeries(psd, (INT4)(f0/delta_f), shrink_length);

  return(psd);
}

/* return a unity response function, for use with calibrated data */
COMPLEX8FrequencySeries *unity_response(LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length)
{
  /* variables */
  COMPLEX8FrequencySeries *response;
  int i;

  /* allocate memory */
  response = XLALCreateCOMPLEX8FrequencySeries("response", &epoch, f0, \
      delta_f, &units, length);

  /* get unity response function */
  for (i = 0; i < length; i++)
  {
    response->data->data[i].re = 1;
    response->data->data[i].im = 0;
  }

  return(response);
}

/* generate response function for LIGO data */
COMPLEX8FrequencySeries *ligo_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length,
    INT4 offset)
{
  /* variables */
  COMPLEX8FrequencySeries *response;
  FrCache *cache = NULL;
  CalibrationUpdateParams calib_params;

  /* apply offset to epoch */
  epoch.gpsSeconds += offset;

  /* allocate memory */
  memset(&calib_params, 0, sizeof(CalibrationUpdateParams));
  response = XLALCreateCOMPLEX8FrequencySeries("response", &epoch, 0, \
      delta_f, &units, length + (INT4)(f0/delta_f));

  /* set ifo */
  calib_params.ifo = ifo;

  /* open calibration frame cache */
  cache = XLALFrImportCache(cache_file);

  /* generate response function */
  LAL_CALL(LALExtractFrameResponse(status, response, cache, &calib_params), \
      status);

  /* destory calibration frame cache */
  XLALFrDestroyCache(cache);

  /* reduce to required band */
  XLALShrinkCOMPLEX8FrequencySeries(response, f0/delta_f, length);

  return(response);
}

/* wrapper to unity_response() and ligo_response() for generating the
 * appropriate response for the given detector */
COMPLEX8FrequencySeries *generate_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length,
    INT4 offset)
{
  /* variables */
  COMPLEX8FrequencySeries *response = NULL;

  if (strncmp(ifo, "G1", 2) == 0)
  {
    /* generate response for GEO */
    response = unity_response(epoch, f0, delta_f, units, length);
  }
  else
  {
    /* generate response function for LIGO */
    response = ligo_response(status, ifo, cache_file, epoch, f0, delta_f, \
        units, length, offset);
  }

  return(response);
}

/* return the frequency mask */
REAL4FrequencySeries *frequency_mask(REAL8 f0,
    REAL8 delta_f,
    INT4 length,
    INT4 bins)
{
  /* counters */
  INT4 i, j;

  /* variables */
  REAL4FrequencySeries *mask;
  LIGOTimeGPS epoch;
  INT4 nBins;
  INT4 numFMin;
  INT4 full_length;

  /* initialise time */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* extra bins */
  nBins = (bins - 1) / 2;
  numFMin = (INT4)(f0 / delta_f);

  /* get length for full frequency band */
  full_length = length + numFMin;

  /* allocate memory for frequency mask */
  mask = XLALCreateREAL4FrequencySeries("mask", &epoch, 0, delta_f, \
      &lalDimensionlessUnit, full_length);

  /* set all values to 1 */
  for (i = 0; i < full_length; i++)
    mask->data->data[i] = 1;

  /* remove multiples of 16 Hz */
  for (i = 0; i < full_length; i += (INT4)(16 / delta_f))
  {
    mask->data->data[i]= 0;

    for (j = 0; j < nBins; j++)
    {
      if ((i + 1 + j) < full_length)
        mask->data->data[i + 1 + j]= 0;
      if ((i - 1 - j) > 0 )
        mask->data->data[i - 1 - j]= 0;
    }
  }

  /* remove multiples of 60 Hz */
  for (i = 0; i < full_length; i += (INT4)(60 / delta_f))
  {
    mask->data->data[i] = 0;

    for (j = 0; j < nBins; j++)
    {
      if ((i + 1 + j) < full_length)
        mask->data->data[i + 1 + j]= 0;
      if ((i - 1 - j) > 0 )
        mask->data->data[i - 1 - j]= 0;
    }
  }

  /* get appropriate band */
  XLALShrinkREAL4FrequencySeries(mask, numFMin, length);

  return(mask);
}

/* zero pad and fft */
COMPLEX8FrequencySeries *zero_pad_and_fft(LALStatus *status,
    REAL4TimeSeries *series,
    REAL8 delta_f,
    INT4 length,
    REAL4Window *window)
{
  /* variables */
  COMPLEX8FrequencySeries *zero_pad;
  RealFFTPlan *plan = NULL;
  SZeroPadAndFFTParameters zero_pad_params;

  /* create fft plan */
  plan = XLALCreateForwardREAL4FFTPlan(2 * series->data->length, 0);

  /* allocate memory */
  zero_pad = XLALCreateCOMPLEX8FrequencySeries("zero_pad", &series->epoch, \
      0, delta_f, &lalDimensionlessUnit, length);

  /* set zeropad parameters */
  zero_pad_params.fftPlan = plan;
  zero_pad_params.window = window->data;
  zero_pad_params.length = 2 * series->data->length;

  /* zero pad and fft */
  LAL_CALL(LALSZeroPadAndFFT(status, zero_pad, series, &zero_pad_params), \
      status);

  /* destroy fft plan */
  XLALDestroyREAL4FFTPlan(plan);

  return(zero_pad);
}

/* generate the cross correlation spectra */
COMPLEX8FrequencySeries *cc_spectrum(LALStatus *status,
    COMPLEX8FrequencySeries *zero_pad_one,
    COMPLEX8FrequencySeries *zero_pad_two,
    COMPLEX8FrequencySeries *response_one,
    COMPLEX8FrequencySeries *response_two,
    REAL4FrequencySeries *opt_filter)
{
  /* variables */
  COMPLEX8FrequencySeries *series;
  StochasticCrossCorrelationCalInput cc_input;

  /* allocate memory */
  series = XLALCreateCOMPLEX8FrequencySeries("cc_spectra", \
      &opt_filter->epoch, opt_filter->f0, opt_filter->deltaF, \
      &lalDimensionlessUnit, opt_filter->data->length);

  /* set inputs */
  cc_input.hBarTildeOne = zero_pad_one;
  cc_input.hBarTildeTwo = zero_pad_two;
  cc_input.responseFunctionOne = response_one;
  cc_input.responseFunctionTwo = response_two;
  cc_input.optimalFilter = opt_filter;

  /* calculate spectrum */
  LAL_CALL(LALStochasticCrossCorrelationSpectrumCal(status, series, \
        &cc_input, 1), status);

  return(series);
}

/* generate cross correlation statistic from cross correlation spectra */
REAL8 cc_statistic(COMPLEX8FrequencySeries *cc_spectra)
{
  /* variables */
  REAL8 cc_stat = 0;
  UINT4 i;

  /* sum up frequencies */
  for (i = 0; i < cc_spectra->data->length; i++)
  {
    cc_stat += cc_spectra->data->data[i].re;
  }

  /* normalise */
  cc_stat *= 2 * cc_spectra->deltaF;

  return(cc_stat);
}

/* main stochastic search */
StochasticTable *stochastic_search(LALStatus *status,
    REAL4TimeSeries *series_one,
    REAL4TimeSeries *series_two,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4Window *window,
    INT4 num_segments,
    INT4 filter_length,
    INT4 segment_duration,
    INT4 segs_in_interval,
    INT4 segment_length,
    INT4 psd_window_duration,
    CHAR *ifo_one,
    CHAR *ifo_two,
    CHAR *channel_one,
    CHAR *channel_two,
    CHAR *calibration_cache_one,
    CHAR *calibration_cache_two,
    INT4 calibration_offset,
    REAL8 f_min,
    REAL8 f_max,
    REAL8 f_ref)
{
  /* counters */
  INT4 i, j, k;

  /* data structures */
  REAL4Vector *cal_psd_one;
  REAL4Vector *cal_psd_two;
  REAL4TimeSeries *segment_one = NULL;
  REAL4TimeSeries *segment_two = NULL;
  COMPLEX8FrequencySeries *response_one = NULL;
  COMPLEX8FrequencySeries *response_two = NULL;
  REAL4FrequencySeries *psd_one = NULL;
  REAL4FrequencySeries *psd_two = NULL;
  REAL4FrequencySeries *inv_psd_one = NULL;
  REAL4FrequencySeries *inv_psd_two = NULL;
  REAL4FrequencySeries *opt_filter = NULL;
  COMPLEX8FrequencySeries *zero_pad_one = NULL;
  COMPLEX8FrequencySeries *zero_pad_two = NULL;
  COMPLEX8FrequencySeries *cc_spectra = NULL;
  LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* variables */
  INT4 num_intervals;
  INT4 segment_shift;
  INT4 middle_segment;
  LIGOTimeGPS seg_epoch;
  LIGOTimeGPS analysis_epoch;
  REAL8 delta_f;

  /* results */
  REAL8 sigma;
  REAL8 y;
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;

  /* debug file */
  CHAR debug_filename[FILENAME_MAX];

  /* initialise analysis_epoch */
  analysis_epoch.gpsSeconds = 0;
  analysis_epoch.gpsNanoSeconds = 0;

  /* calculate number of intervals, and required shift to get to next
   * interval */
  if (overlap_hann_flag)
  {
    num_intervals = (2 * (num_segments - 2)) - 1;
    segment_shift = segment_duration / 2;
  }
  else
  {
    num_intervals = num_segments - 2;
    segment_shift = segment_duration;
  }

  /* get middle segment number */
  middle_segment = (segs_in_interval - 1) / 2;

  /* get delta_f for optimal filter */
  delta_f = 1./(REAL8)psd_window_duration;

  /* allocate memory for calibrated PSDs */
  cal_psd_one = XLALCreateREAL4Vector(filter_length);
  cal_psd_two = XLALCreateREAL4Vector(filter_length);

  if (vrbflg)
    fprintf(stdout, "Starting analysis loop...\n");

  /* loop over intervals */
  for (j = 0; j < num_intervals; j++)
  {	
    /* initialize average PSDs */
    for (i = 0; i < filter_length; i++)
    {
      cal_psd_one->data[i] = 0;
      cal_psd_two->data[i] = 0;
    }

    /* loop over segments in the interval */
    for (k = 0; k < segs_in_interval; k++)
    {
      /* set segment start time */
      seg_epoch = series_one->epoch;
      XLALGPSAdd(&seg_epoch, j * segment_shift + k * segment_duration);

      /* is this the analysis segment? */
      if (k == middle_segment)
        analysis_epoch = seg_epoch;

      if (vrbflg)
      {
        fprintf(stdout, "Request data at GPS time %d\n", \
            seg_epoch.gpsSeconds);
      }

      /* cut segments from series */
      segment_one = cut_time_series(series_one, seg_epoch, segment_duration);
      segment_two = cut_time_series(series_two, seg_epoch, segment_duration);

      /* save intermediate products */
      if (debug_flag)
      {
        snprintf(debug_filename, FILENAME_MAX, "segment_1-%d.dat", \
            seg_epoch.gpsSeconds);
        LALSPrintTimeSeries(segment_one, debug_filename);
        snprintf(debug_filename, FILENAME_MAX, "segment_2-%d.dat", \
            seg_epoch.gpsSeconds);
        LALSPrintTimeSeries(segment_two, debug_filename);
      }

      /* compute response */
      response_one = generate_response(status, ifo_one, calibration_cache_one, \
          seg_epoch, f_min, delta_f, countPerAttoStrain, filter_length, \
          calibration_offset);
      response_two = generate_response(status, ifo_two, calibration_cache_two, \
          seg_epoch, f_min, delta_f, countPerAttoStrain, filter_length, \
          calibration_offset);

      /* save intermediate products */
      if (debug_flag)
      {
        snprintf(debug_filename, FILENAME_MAX, "response_1-%d.dat", \
            seg_epoch.gpsSeconds + calibration_offset);
        LALCPrintFrequencySeries(response_one, debug_filename);
        snprintf(debug_filename, FILENAME_MAX, "response_2-%d.dat", \
            seg_epoch.gpsSeconds + calibration_offset);
        LALCPrintFrequencySeries(response_two, debug_filename);
      }

      /* check if on middle segment and if we want to include this in
       * the analysis */
      if ((k == middle_segment) && (!middle_segment_flag))
      {
        if (vrbflg)
          fprintf(stdout, "Ignoring middle segment..\n");
      }
      else
      {
        if (vrbflg)
          fprintf(stdout, "Estimating PSDs...\n");

        /* compute uncalibrated PSDs */
        psd_one = estimate_psd(segment_one, f_min, filter_length, \
            psd_window_duration);
        psd_two = estimate_psd(segment_two, f_min, filter_length, \
            psd_window_duration);

        if (vrbflg)
          fprintf(stdout, "Generating inverse noise...\n");

        /* compute inverse calibrate PSDs */
        inv_psd_one = inverse_noise(status, psd_one, response_one);
        inv_psd_two = inverse_noise(status, psd_two, response_two);

        /* sum over calibrated PSDs for average */
        for (i = 0; i < filter_length; i++)
        {
          cal_psd_one->data[i] += 1. / inv_psd_one->data->data[i];
          cal_psd_two->data[i] += 1. / inv_psd_two->data->data[i];
        }
      }
    }

    /* average calibrated PSDs and take inverse */
    for (i = 0; i < filter_length; i++)
    {
      /* average */
      if (!middle_segment_flag)
      {
        cal_psd_one->data[i] /= (REAL4)(segs_in_interval - 1);
        cal_psd_two->data[i] /= (REAL4)(segs_in_interval - 1);
      }
      else
      {
        cal_psd_one->data[i] /= (REAL4)segs_in_interval;
        cal_psd_two->data[i] /= (REAL4)segs_in_interval;
      }
      /* take inverse */
      inv_psd_one->data->data[i] = 1. / cal_psd_one->data[i];
      inv_psd_two->data->data[i] = 1. / cal_psd_two->data[i];
    }

    /* save intermediate products */
    if (debug_flag)
    {
      snprintf(debug_filename, FILENAME_MAX, "inv_psd_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(inv_psd_one, debug_filename);
      snprintf(debug_filename, FILENAME_MAX, "inv_psd_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(inv_psd_two, debug_filename);
    }

    if (vrbflg)
      fprintf(stdout, "Generating optimal filter...\n");

    /* build optimal filter */
    opt_filter = optimal_filter(status, overlap, omega, inv_psd_one, \
        inv_psd_two, window, &sigma, f_ref, segment_duration);

    /* save intermediate products */
    if (debug_flag)
    {
      snprintf(debug_filename, FILENAME_MAX, "opt_filter-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(opt_filter, debug_filename);
    }          

    if (vrbflg)
    {
      fprintf(stdout, "Analysing segment at GPS %d\n", \
          analysis_epoch.gpsSeconds);
    }

    /* cut analysis segment from full series */
    segment_one = cut_time_series(series_one, analysis_epoch, \
        segment_duration);
    segment_two = cut_time_series(series_two, analysis_epoch, \
        segment_duration);

    /* save intermediate products */
    if (debug_flag)
    {
      snprintf(debug_filename, FILENAME_MAX, "analysis_segment_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintTimeSeries(segment_one, debug_filename);
      snprintf(debug_filename, FILENAME_MAX, "analysis_segment_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintTimeSeries(segment_two, debug_filename);
    }

    if (vrbflg)
      fprintf(stdout, "Performing zero pad and FFT...\n");

    /* zero pad and fft */
    zero_pad_one = zero_pad_and_fft(status, segment_one, delta_f, \
        segment_length + 1, window);
    zero_pad_two = zero_pad_and_fft(status, segment_two, delta_f, \
        segment_length + 1, window);

    /* save intermediate products */
    if (debug_flag)
    {
      snprintf(debug_filename, FILENAME_MAX, "zero_pad_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALCPrintFrequencySeries(zero_pad_one, debug_filename);
      snprintf(debug_filename, FILENAME_MAX, "zero_pad_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALCPrintFrequencySeries(zero_pad_two, debug_filename);
    }    

    if (vrbflg)
      fprintf(stdout, "Calculating cross correlation spectrum...\n");

    /* get response functions for analysis segments */
    response_one = generate_response(status, ifo_one, calibration_cache_one, \
        analysis_epoch, f_min, delta_f, countPerAttoStrain, filter_length, \
        calibration_offset);
    response_two = generate_response(status, ifo_two, calibration_cache_two, \
        analysis_epoch, f_min, delta_f, countPerAttoStrain, filter_length, \
        calibration_offset);

    /* calculate cc spectrum */
    cc_spectra = cc_spectrum(status, zero_pad_one, zero_pad_two, \
        response_one, response_two, opt_filter);

    if (cc_spectra_flag)
    {
      /* save out cc spectra as frame */
      if (vrbflg)
        fprintf(stdout, "Saving ccSpectra to frame...\n");
      write_ccspectra_frame(cc_spectra, ifo_one, ifo_two, \
          analysis_epoch, segment_duration);
    }

    /* cc statistic */
    y = cc_statistic(cc_spectra);

    /* display results */
    if (vrbflg)
    {
      fprintf(stdout, "Interval %d\n", j + 1);
      fprintf(stdout, "  GPS time  = %d\n", analysis_epoch.gpsSeconds);
      fprintf(stdout, "  y         = %e\n", y);
      fprintf(stdout, "  sigma     = %e\n", sigma);
    }

    /* allocate memory for table */
    if (!stochHead)
    {
      stochHead = thisStoch = (StochasticTable *) \
                  LALCalloc(1, sizeof(StochasticTable));
    }
    else
    {
      thisStoch = thisStoch->next = (StochasticTable *) \
                  LALCalloc(1, sizeof(StochasticTable));
    }

    /* populate columns */
    snprintf(thisStoch->ifo_one, LIGOMETA_IFO_MAX, ifo_one);
    snprintf(thisStoch->ifo_two, LIGOMETA_IFO_MAX, ifo_two);
    snprintf(thisStoch->channel_one, LIGOMETA_CHANNEL_MAX, channel_one);
    snprintf(thisStoch->channel_two, LIGOMETA_CHANNEL_MAX, channel_two);
    thisStoch->start_time.gpsSeconds = analysis_epoch.gpsSeconds;
    thisStoch->start_time.gpsNanoSeconds = analysis_epoch.gpsNanoSeconds;
    thisStoch->duration.gpsSeconds = segment_duration;
    thisStoch->duration.gpsNanoSeconds = 0;
    thisStoch->f_min = f_min;
    thisStoch->f_max = f_max;
    thisStoch->cc_stat = y;
    thisStoch->cc_sigma = sigma;
  }

  /* clean up */
  XLALDestroyREAL4Vector(cal_psd_one);
  XLALDestroyREAL4Vector(cal_psd_two);
  XLALDestroyREAL4TimeSeries(segment_one);
  XLALDestroyREAL4TimeSeries(segment_two);
  XLALDestroyCOMPLEX8FrequencySeries(response_one);
  XLALDestroyCOMPLEX8FrequencySeries(response_two);
  XLALDestroyREAL4FrequencySeries(psd_one);
  XLALDestroyREAL4FrequencySeries(psd_two);
  XLALDestroyREAL4FrequencySeries(inv_psd_one);
  XLALDestroyREAL4FrequencySeries(inv_psd_two);
  XLALDestroyREAL4FrequencySeries(opt_filter);
  XLALDestroyCOMPLEX8FrequencySeries(zero_pad_one);
  XLALDestroyCOMPLEX8FrequencySeries(zero_pad_two);
  XLALDestroyCOMPLEX8FrequencySeries(cc_spectra);

  return(stochHead);
}

/*
 * vim: et
 */
