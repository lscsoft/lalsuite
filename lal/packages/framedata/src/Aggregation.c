/*
 * Aggregation.c - online frame data aggregation routines
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2009 Adam Mercer
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Aggregation.h>
#include <lal/XLALError.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/TimeSeries.h>


/* return frame gps start time for given gps time */
LIGOTimeGPS *XLALAggregationFrameStart(LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationFrameStart";

  /* declare variables */
  LIGOTimeGPS *start;

  /* check arguments */
  if (!gps)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* allocate memory */
  start = XLALCalloc(1, sizeof(*start));
  if (start == NULL)
  {
    /* failed to allocate memory for start */
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }

  /* determine frame start time, multiple of ONLINE_FRAME_DURATION */
  start->gpsSeconds = (INT4)floor(gps->gpsSeconds / ONLINE_FRAME_DURATION) * \
                     ONLINE_FRAME_DURATION;
  start->gpsNanoSeconds = 0;

  return start;
}


/* return frame type for given ifo */
CHAR *XLALAggregationFrameType(CHAR *ifo)
{
  static const char *func = "XLALAggregationFrameType";

  /* declare variables */
  static CHAR type[LIGOMETA_TYPE_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine type */
  if (strncmp(ifo, "G1", LIGOMETA_IFO_MAX) == 0)
  {
    /* FIXME geo - currently undefined */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }
  else if (strncmp(ifo, "V1", LIGOMETA_IFO_MAX) == 0)
  {
    /* FIXME virgo - currently undefined */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }
  else if ((strncmp(ifo, "H1", LIGOMETA_IFO_MAX) == 0) || \
      (strncmp(ifo, "H2", LIGOMETA_IFO_MAX) == 0) || \
      (strncmp(ifo, "L1", LIGOMETA_IFO_MAX) == 0))
  {
    /* ligo */
    LALSnprintf(type, LIGOMETA_TYPE_MAX, "%s_DMT_C00_L2", ifo);
  }
  else
  {
    /* unsupported ifo */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  return type;
}


/* return path for given ifo and gps time */
CHAR *XLALAggregationDirectoryPath(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationDirectoryPath";

  /* declare variables */
  CHAR *base_dir;
  CHAR *type;
  LIGOTimeGPS *frame_start;
  INT4 gps_dir;
  static CHAR directory[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* get base directory from environment */
  base_dir = getenv("ONLINEHOFT");
  if (base_dir == NULL)
  {
    /* ONLINEHOFT environment variable not set */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* determine gps directory name, frame start must be multiple of
   * ONLINE_FRAME_DURATION */
  frame_start = XLALAggregationFrameStart(gps);
  gps_dir = (INT4)floor(frame_start->gpsSeconds / 100000);

  /* construct directory */
  LALSnprintf(directory, FILENAME_MAX, "%s/%c-%s-%d", base_dir, ifo[0], \
      type, gps_dir);

  /* free memory */
  XLALFree(frame_start);

  return directory;
}


/* return frame filename for given ifo and gps time */
CHAR *XLALAggregationFrameFilename(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationFrameFilename";

  /* declare variables */
  CHAR *type;
  LIGOTimeGPS *frame_start;
  static CHAR filename[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* determine gps start time for frame */
  frame_start = XLALAggregationFrameStart(gps);

  /* construct frame filename */
  LALSnprintf(filename, FILENAME_MAX, "%c-%s-%d-%d.gwf", ifo[0], type, \
      frame_start->gpsSeconds, ONLINE_FRAME_DURATION);

  /* free memory */
  XLALFree(frame_start);

  return filename;
}


/* return full path to frame for given ifo and gps time */
CHAR *XLALAggregationFramePathFilename(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationFramePathFilename";

  /* declare variables */
  CHAR *directory;
  CHAR *filename;
  static CHAR path[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine directory */
  directory = XLALAggregationDirectoryPath(ifo, gps);
  if (directory == NULL)
  {
    /* failed to determine directory */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* determine frame filename */
  filename = XLALAggregationFrameFilename(ifo, gps);
  if (filename == NULL)
  {
    /* failed to determine filename */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* construct path */
  LALSnprintf(path, FILENAME_MAX, "%s/%s", directory, filename);

  return path;
}


/* return url to frame for a given ifo and gps time */
CHAR *XLALAggregationFrameURL(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationFrameURL";

  /* declare variables */
  CHAR *path;
  static CHAR url[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine path */
  path = XLALAggregationFramePathFilename(ifo, gps);
  if (path == NULL)
  {
    /* failed to determine path */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* construct url */
  LALSnprintf(url, FILENAME_MAX, "file://localhost%s", path);

  return url;
}


/* return frame cache given ifo, gps time, and length */
FrCache *XLALAggregationFrameCache(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length)
{
  static const char *func = "XLALAggregationFrameCache";

  /* declare variables */
  LIGOTimeGPS gps;
  LIGOTimeGPS *frame_start;
  LIGOTimeGPS *last_frame_start;
  CHAR *type;
  CHAR *url;
  FrCache *cache;
  int i = 0;
  INT4 frame_duration;
  INT4 num_frames;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine number of frames */
  gps.gpsSeconds = start->gpsSeconds + length;
  gps.gpsNanoSeconds = start->gpsNanoSeconds;
  frame_start = XLALAggregationFrameStart(start);
  last_frame_start = XLALAggregationFrameStart(&gps);
  frame_duration = (last_frame_start->gpsSeconds + ONLINE_FRAME_DURATION) - \
                   frame_start->gpsSeconds;
  num_frames = frame_duration / ONLINE_FRAME_DURATION;

  /* free memory */
  XLALFree(frame_start);
  XLALFree(last_frame_start);

  /* initilise cache */
  cache = XLALCalloc(1, sizeof(*cache));
  if (cache == NULL)
  {
    /* failed to allocate memory for cache */
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }
  cache->numFrameFiles = num_frames;
  cache->frameFiles = XLALCalloc(num_frames, sizeof(*cache->frameFiles));
  if (cache->frameFiles == NULL)
  {
    /* failed to allocate memory for cache->frameFiles */
    XLALFree(cache);
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLALFrDestroyCache(cache);
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* initialise gps */
  gps.gpsSeconds = start->gpsSeconds;
  gps.gpsNanoSeconds = start->gpsNanoSeconds;

  /* determine list of frames */
  for (i = 0; i < num_frames; i++)
  {
    /* declare variables */
    FrStat *file = cache->frameFiles + i;

    /* increment gps */
    gps.gpsSeconds = start->gpsSeconds + (i * ONLINE_FRAME_DURATION);
    gps.gpsNanoSeconds = start->gpsNanoSeconds;

    /* determine url */
    url = XLALAggregationFrameURL(ifo, &gps);
    if (url == NULL)
    {
      /* failed to determine url */
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_EINVAL);
    }

    /* determine frame start */
    frame_start = XLALAggregationFrameStart(&gps);

    /* allocate memory for cache entry */
    file->source = XLALMalloc(strlen(ifo) + 1);
    if (!file->source)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }
    file->description = XLALMalloc(strlen(type) + 1);
    if (!file->description)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }
    file->url = XLALMalloc(strlen(url) + 1);
    if (!file->url)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }

    /* add frame to cache */
    strcpy(file->source, ifo);
    strcpy(file->description, type);
    file->startTime = frame_start->gpsSeconds;
    file->duration = ONLINE_FRAME_DURATION;
    strcpy(file->url, url);

    /* free memory */
    XLALFree(frame_start);
  }

  return cache;
}


/* return frame stream for given ifo, gps time, and length */
FrStream *XLALAggregationFrameStream(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length)
{
  static const char *func = "XLALAggregationFrameStream";

  /* declare variables */
  FrCache *cache;
  FrStream *stream;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* get frame cache */
  cache = XLALAggregationFrameCache(ifo, start, length);
  if (cache == NULL)
  {
    /* failed to get cache */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* open cache as stream */
  stream = XLALFrCacheOpen(cache);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLALFrDestroyCache(cache);
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* destroy cache */
  XLALFrDestroyCache(cache);

  return stream;
}


/* return strain data time series for given ifo, gps time, and length */
REAL8TimeSeries *XLALAggregationStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length)
{
  static const char *func = "XLALAggregationStrainData";

  /* declare variables */
  FrStream *stream;
  REAL8TimeSeries *series;
  size_t num_points;
  CHAR channel[LIGOMETA_CHANNEL_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine number of data points */
  num_points = length * ONLINE_SAMPLE_RATE;

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, length);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* get strain data time series */
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "%s:%s", ifo, \
      ONLINE_STRAIN_CHANNEL);
  series = XLALFrReadREAL8TimeSeries(stream, channel, start, \
      (REAL8)length, num_points);
  if (series == NULL)
  {
    /* failed to read data */
    XLALFrClose(stream);
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* close stream */
  XLALFrClose(stream);

  return series;
}


/* return data quality vector time series for given ifo, gps time and
 * length */
INT4TimeSeries *XLALAggregationDQVector(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length)
{
  static const char *func = "XLALAggregationDQVector";

  /* declare variables */
  FrStream *stream;
  INT4TimeSeries *series;
  size_t num_points;
  CHAR channel[LIGOMETA_CHANNEL_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine number of data points */
  num_points = length * ONLINE_SAMPLE_RATE;

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, length);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* get data quality vector time series */
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "%s:%s", ifo, \
      ONLINE_DQ_VECTOR);
  series = XLALFrReadINT4TimeSeries(stream, channel, start, \
      (REAL8)length, num_points);
  if (series == NULL)
  {
    /* failed to read data */
    XLALFrClose(stream);
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* close stream */
  XLALFrClose(stream);

  return series;
}


/* return state vector time series for given ifo, gps time, and length */
INT4TimeSeries *XLALAggregationStateVector(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length)
{
  static const char *func = "XLALAggregationStateVector";

  /* declare variables */
  FrStream *stream;
  REAL4TimeSeries *state;
  INT4TimeSeries *series;
  size_t num_points;
  CHAR channel[LIGOMETA_CHANNEL_MAX];
  UINT4 i;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* determine number of data points */
  num_points = length * ONLINE_SAMPLE_RATE;

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, length);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* get state vector time series */
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "%s:%s", ifo, \
      ONLINE_STATE_VECTOR);
  state = XLALFrReadREAL4TimeSeries(stream, channel, start, \
      (REAL8)length, num_points);
  if (state == NULL)
  {
    /* failed to read data */
    XLALFrClose(stream);
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* initialise series */
  series = XLALCreateINT4TimeSeries(state->name, &state->epoch, \
      state->f0, state->deltaT, &state->sampleUnits, \
      state->data->length);

  /* cast state vector to INT4 */
  for (i = 0; i < state->data->length; i++)
    series->data->data[i] = (INT4)state->data->data[i];

  /* destroy state vector */
  XLALDestroyREAL4TimeSeries(state);

  /* close stream */
  XLALFrClose(stream);

  return series;
}


/* return strain data time series for given ifo, gps time, length, and
 * data quality vector bitmask */
REAL8TimeSeries *XLALAggregationDQStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length,
    INT4 dq_bitmask)
{
  static const char *func = "XLALAggregationDQStrainData";

  /* declare variables */
  INT4TimeSeries *dq_vector;
  REAL8TimeSeries *series;
  UINT4 i;

  /* checkout arguments */
  if (!ifo)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* get data quality vector */
  dq_vector = XLALAggregationDQVector(ifo, start, length);
  if (dq_vector == NULL)
  {
    /* failed to get data quality vector */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* check for required bitmask */
  for (i = 0; i < dq_vector->data->length; i++)
  {
    if ((dq_vector->data->data[i] & dq_bitmask) != dq_bitmask)
    {
      /* invalid bitmask */
      XLALDestroyINT4TimeSeries(dq_vector);
      XLAL_ERROR_NULL(func, XLAL_EINVAL);
    }
  }

  /* destroy data quality vector */
  XLALDestroyINT4TimeSeries(dq_vector);

  /* get strain data time series */
  series = XLALAggregationStrainData(ifo, start, length);
  if (series == NULL)
  {
    /* failed to get strain data time series */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  return series;
}
