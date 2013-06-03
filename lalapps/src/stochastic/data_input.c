/*
 * data_input.c - SGWB Standalone Analysis Pipeline
 *              - Data Input Functions
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

#include "data_input.h"


/* externally declared variables */
extern int high_pass_flag;
extern int vrbflg;

/* read a LIGO time series */
REAL4TimeSeries *get_ligo_data(FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  REAL4TimeSeries *series;
  size_t length;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  series = XLALCreateREAL4TimeSeries(channel, &start, 0, 0, \
      &lalADCCountUnit, 0);

  if (vrbflg)
    fprintf(stdout, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  XLALFrGetREAL4TimeSeriesMetadata(series, stream);

  if (vrbflg)
    fprintf(stdout, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = floor((XLALGPSDiff(&end, &start) / series->deltaT) + 0.5);
  XLALResizeREAL4TimeSeries(series, 0, length);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  XLALFrSeek(stream, &start);
  XLALFrGetREAL4TimeSeries(series, stream);

  return(series);
}

/* read and high pass filter a GEO time series */
REAL4TimeSeries *get_geo_data(FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 hpf_order,
    REAL8 hpf_frequency,
    REAL8 hpf_attenuation)
{
  /* variables */
  PassBandParamStruc high_pass_params;
  REAL4TimeSeries *series;
  REAL8TimeSeries *geo;
  size_t length;
  size_t i;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  geo = XLALCreateREAL8TimeSeries(channel, &start, 0, 0, \
        &lalADCCountUnit, 0);

  if (vrbflg)
    fprintf(stdout, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  XLALFrGetREAL8TimeSeriesMetadata(geo, stream);

  if (vrbflg)
    fprintf(stdout, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = floor((XLALGPSDiff(&end, &start) / geo->deltaT) + 0.5);
  XLALResizeREAL8TimeSeries(geo, 0, length);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  XLALFrSeek(stream, &start);
  XLALFrGetREAL8TimeSeries(geo, stream);

  if (vrbflg)
    fprintf(stdout, "High pass filtering \"%s\"...\n", channel);

  /* high pass filter before casting to a REAL4 */
  high_pass_params.nMax = hpf_order;
  high_pass_params.f1 = -1;
  high_pass_params.f2 = hpf_frequency;
  high_pass_params.a1 = -1;
  high_pass_params.a2 = hpf_attenuation;
  XLALButterworthREAL8TimeSeries(geo, &high_pass_params);

  if (vrbflg)
    fprintf(stdout, "Casting \"%s\" as a REAL4...\n", channel);

  /* cast as a REAL4 */
  series = XLALCreateREAL4TimeSeries(geo->name, &geo->epoch, geo->f0, \
      geo->deltaT, &geo->sampleUnits, geo->data->length);
  for (i = 0; i < series->data->length; i++)
    series->data->data[i] = (REAL4)geo->data->data[i];

  /* destroy geo series */
  XLALDestroyREAL8TimeSeries(geo);

  return(series);
}

/* read a time series */
REAL4TimeSeries *get_time_series(CHAR *ifo,
    CHAR *cache_file,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 resample_rate,
    INT4 hpf_order,
    REAL8 hpf_frequency,
    REAL8 hpf_attenuation,
    INT4 geo_hpf_order,
    REAL8 geo_hpf_frequency,
    REAL8 geo_hpf_attenuation,
    INT4 buffer)
{
  /* variables */
  REAL4TimeSeries *series;
  FrStream *stream = NULL;
  LALCache *cache = NULL;
  size_t length;
  PassBandParamStruc high_pass_params;
  int mode = LAL_FR_VERBOSE_MODE;

  /* apply resample buffer if required */
  if (buffer)
  {
    start.gpsSeconds -= buffer;
    end.gpsSeconds += buffer;
  }

  if (vrbflg)
    fprintf(stdout, "Opening frame cache \"%s\"...\n", cache_file);

  /* open frame stream */
  cache = XLALCacheImport(cache_file);
  stream = XLALFrCacheOpen(cache);
  XLALDestroyCache(cache);

  /* turn on checking for missing data */
  XLALFrSetMode(stream, mode);

  /* get the data */
  if (strncmp(ifo, "G1", 2) == 0)
  {
    series = get_geo_data(stream, channel, start, end, geo_hpf_order, \
        geo_hpf_frequency, geo_hpf_attenuation);
  }
  else
    series = get_ligo_data(stream, channel, start, end);

  /* check for missing data */
  if (stream->state & LAL_FR_GAP)
  {
    fprintf(stderr, "Gap in data detected between GPS times %d s and %d s\n", \
        start.gpsSeconds, end.gpsSeconds);
    XLALDestroyREAL4TimeSeries(series);
    exit(1);
  }

  /* clean up */
  XLALFrClose(stream);

  /* resample if required */
  if (resample_rate)
  {
    if (vrbflg)
      fprintf(stdout, "Resampling to %d Hz...\n", resample_rate);

    /* resample */
    XLALResampleREAL4TimeSeries(series, 1./resample_rate);
  }

  /* high pass fitering */
  if (high_pass_flag)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Applying high pass filter to \"%s\"...\n", \
          series->name);
    }

    /* set high pass filter parameters */
    high_pass_params.nMax = hpf_order;
    high_pass_params.f1 = -1;
    high_pass_params.f2 = hpf_frequency;
    high_pass_params.a1 = -1;
    high_pass_params.a2 = hpf_attenuation;

    /* high pass filter */
    XLALButterworthREAL4TimeSeries(series, &high_pass_params);
  }

  /* remove resample buffer */
  if (buffer)
  {
    /* recover original start and end times */
    start.gpsSeconds += buffer;
    end.gpsSeconds -= buffer;

    /* calculate required length */
    length = floor((XLALGPSDiff(&end, &start) / series->deltaT) + 0.5);

    /* remove resample buffer */
    XLALShrinkREAL4TimeSeries(series, buffer / series->deltaT, length);
  }

  return(series);
}

/*
 * vim: et
 */
