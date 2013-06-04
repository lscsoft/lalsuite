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
#include <unistd.h>
#include <sys/stat.h>

#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/Aggregation.h>
#include <lal/XLALError.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/TimeSeries.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>


/* return frame gps start time for given gps time */
LIGOTimeGPS *XLALAggregationFrameStart(LIGOTimeGPS *gps)
{
  /* declare variables */
  LIGOTimeGPS *start;

  /* check arguments */
  if (!gps)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* allocate memory */
  start = XLALCalloc(1, sizeof(*start));
  if (start == NULL)
  {
    /* failed to allocate memory for start */
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  /* determine frame start time, multiple of LAL_ONLINE_FRAME_DURATION */
  start->gpsSeconds = (INT4)floor(gps->gpsSeconds / LAL_ONLINE_FRAME_DURATION) * \
                      LAL_ONLINE_FRAME_DURATION;
  start->gpsNanoSeconds = 0;

  return start;
}


/* return frame type for given ifo */
CHAR *XLALAggregationFrameType(CHAR *ifo)
{
  /* declare variables */
  static CHAR type[LAL_AGGREGATION_TYPE_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine type */
  if (strncmp(ifo, "G1", LAL_AGGREGATION_IFO_MAX) == 0)
  {
    /* FIXME geo - currently undefined */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }
  else if (strncmp(ifo, "V1", LAL_AGGREGATION_IFO_MAX) == 0)
  {
    /* virgo */
    snprintf(type, LAL_AGGREGATION_TYPE_MAX, "%s_DMT_HREC", ifo);
  }
  else if ((strncmp(ifo, "H1", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "H2", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "L1", LAL_AGGREGATION_IFO_MAX) == 0))
  {
    /* ligo */
    snprintf(type, LAL_AGGREGATION_TYPE_MAX, "%s_DMT_C00_L2", ifo);
  }
  else
  {
    /* unsupported ifo */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  return type;
}


/* return path for given ifo and gps time */
CHAR *XLALAggregationDirectoryPath(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  /* declare variables */
  CHAR *base_dir;
  CHAR *type;
  LIGOTimeGPS *frame_start;
  INT4 gps_dir;
  static CHAR directory[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get base directory from environment */
  base_dir = getenv("ONLINEHOFT");
  if (base_dir == NULL)
  {
    /* ONLINEHOFT environment variable not set */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* determine gps directory name, frame start must be multiple of
   * LAL_ONLINE_FRAME_DURATION */
  frame_start = XLALAggregationFrameStart(gps);
  gps_dir = (INT4)floor(frame_start->gpsSeconds / 100000);

  /* construct directory */
  snprintf(directory, FILENAME_MAX, "%s/%s/%c-%s-%d", base_dir, ifo, \
      ifo[0], type, gps_dir);

  /* free memory */
  XLALFree(frame_start);

  return directory;
}


/* return frame filename for given ifo and gps time */
CHAR *XLALAggregationFrameFilename(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  /* declare variables */
  CHAR *type;
  LIGOTimeGPS *frame_start;
  static CHAR filename[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* determine gps start time for frame */
  frame_start = XLALAggregationFrameStart(gps);

  /* construct frame filename */
  snprintf(filename, FILENAME_MAX, "%c-%s-%d-%d.gwf", ifo[0], type, \
      frame_start->gpsSeconds, LAL_ONLINE_FRAME_DURATION);

  /* free memory */
  XLALFree(frame_start);

  return filename;
}


/* return full path to frame for given ifo and gps time */
CHAR *XLALAggregationFramePathFilename(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  /* declare variables */
  CHAR *directory;
  CHAR *filename;
  static CHAR path[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine directory */
  directory = XLALAggregationDirectoryPath(ifo, gps);
  if (directory == NULL)
  {
    /* failed to determine directory */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* determine frame filename */
  filename = XLALAggregationFrameFilename(ifo, gps);
  if (filename == NULL)
  {
    /* failed to determine filename */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* construct path */
  snprintf(path, FILENAME_MAX, "%s/%s", directory, filename);

  return path;
}


/* return url to frame for a given ifo and gps time */
CHAR *XLALAggregationFrameURL(CHAR *ifo,
    LIGOTimeGPS *gps)
{
  /* declare variables */
  CHAR *path;
  static CHAR url[FILENAME_MAX];

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!gps)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine path */
  path = XLALAggregationFramePathFilename(ifo, gps);
  if (path == NULL)
  {
    /* failed to determine path */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* construct url */
  snprintf(url, FILENAME_MAX, "file://localhost%s", path);

  return url;
}

/* return gps time of latest frame written to disk */
LIGOTimeGPS *XLALAggregationLatestGPS(CHAR *ifo)
{
  /* declare variables */
  CHAR *base_dir;
  CHAR path[FILENAME_MAX];
  struct stat file_status;
  FILE *file_ptr;
  static LIGOTimeGPS gps = {0, 0};
  int i;

  /* checkout arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine directory from environment */
  base_dir = getenv("ONLINEHOFT");
  if (base_dir == NULL)
  {
    /* ONLINEHOFT environment variable not set */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* determine path to latest file */
  snprintf(path, FILENAME_MAX, "%s/%s/latest", base_dir, ifo);

  /* check file status */
  if (stat(path, &file_status) == -1)
  {
    /* failed to find file */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* open latest file */
  file_ptr = fopen(path, "r");
  if (file_ptr == NULL)
  {
    /* failed to open file */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* determine gps time of latest frame file written */
  i = fscanf(file_ptr, "%d", &gps.gpsSeconds);
  if (i != 1)
  {
    /* failed to get latest gps time */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* close file */
  fclose(file_ptr);

  return &gps;
}


/* return frame cache given ifo, gps time, and duration */
LALCache *XLALAggregationFrameCache(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LIGOTimeGPS gps;
  LIGOTimeGPS *frame_start;
  LIGOTimeGPS *last_frame_start;
  CHAR *type;
  CHAR *url;
  LALCache *cache;
  int i = 0;
  INT4 frame_duration;
  INT4 num_frames;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* determine number of frames */
  gps.gpsSeconds = start->gpsSeconds + (INT4)floor(duration);
  gps.gpsNanoSeconds = start->gpsNanoSeconds + \
                       (INT4)((duration - floor(duration)) * pow(10, 9));
  frame_start = XLALAggregationFrameStart(start);
  last_frame_start = XLALAggregationFrameStart(&gps);
  frame_duration = (last_frame_start->gpsSeconds + LAL_ONLINE_FRAME_DURATION) - \
                   frame_start->gpsSeconds;
  num_frames = frame_duration / LAL_ONLINE_FRAME_DURATION;

  /* free memory */
  XLALFree(frame_start);
  XLALFree(last_frame_start);

  /* initilise cache */
  cache = XLALCalloc(1, sizeof(*cache));
  if (cache == NULL)
  {
    /* failed to allocate memory for cache */
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  cache->length = num_frames;
  cache->list = XLALCalloc(num_frames, sizeof(*cache->list));
  if (cache->list == NULL)
  {
    /* failed to allocate memory for cache->list */
    XLALFree(cache);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  /* determine type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    /* failed to determine type */
    XLALDestroyCache(cache);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* initialise gps */
  gps.gpsSeconds = start->gpsSeconds;
  gps.gpsNanoSeconds = start->gpsNanoSeconds;

  /* determine list of frames */
  for (i = 0; i < num_frames; i++)
  {
    /* declare variables */
    LALCacheEntry *file = cache->list + i;

    /* increment gps */
    gps.gpsSeconds = start->gpsSeconds + (i * LAL_ONLINE_FRAME_DURATION);
    gps.gpsNanoSeconds = start->gpsNanoSeconds;

    /* determine url */
    url = XLALAggregationFrameURL(ifo, &gps);
    if (url == NULL)
    {
      /* failed to determine url */
      XLALDestroyCache(cache);
      XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* determine frame start */
    frame_start = XLALAggregationFrameStart(&gps);

    /* allocate memory for cache entry */
    file->src = XLALMalloc(strlen(ifo) + 1);
    if (!file->src)
    {
      XLALDestroyCache(cache);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    file->dsc = XLALMalloc(strlen(type) + 1);
    if (!file->dsc)
    {
      XLALDestroyCache(cache);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    file->url = XLALMalloc(strlen(url) + 1);
    if (!file->url)
    {
      XLALDestroyCache(cache);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }

    /* add frame to cache */
    strcpy(file->src, ifo);
    strcpy(file->dsc, type);
    file->t0 = frame_start->gpsSeconds;
    file->dt = LAL_ONLINE_FRAME_DURATION;
    strcpy(file->url, url);

    /* free memory */
    XLALFree(frame_start);
  }

  return cache;
}


/* return frame stream for given ifo, gps time, and duration */
LALFrStream *XLALAggregationFrameStream(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LALCache *cache;
  LALFrStream *stream;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get frame cache */
  cache = XLALAggregationFrameCache(ifo, start, duration);
  if (cache == NULL)
  {
    /* failed to get cache */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* open cache as stream */
  stream = XLALFrStreamCacheOpen(cache);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLALDestroyCache(cache);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* destroy cache */
  XLALDestroyCache(cache);

  return stream;
}


/* return strain data time series for given ifo, gps time, and duration */
REAL8TimeSeries *XLALAggregationStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LALFrStream *stream;
  REAL8TimeSeries *series;
  CHAR channel[LAL_AGGREGATION_CHANNEL_MAX];
  LIGOTimeGPS time_now;
  LIGOTimeGPS *latest;
  INT4 end_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get current gps time */
  if (XLALGPSTimeNow(&time_now) == NULL)
  {
    /* failed to get current time */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* check that requested data is not in the future */
  if (XLALGPSCmp(&time_now, start) == -1)
  {
    /* requested time in the future */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* determine gps time of latest frame file written */
  latest = XLALAggregationLatestGPS(ifo);
  if (latest == NULL)
  {
    /* failed to determine gps time */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* get end time of requested data */
  end_time = start->gpsSeconds + (INT4)round(duration);

  /* check that requested data has been written */
  if (latest->gpsSeconds < end_time)
  {
    /* requested data has not been written yet */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get channel name */
  if (strncmp(ifo, "V1", LAL_AGGREGATION_IFO_MAX) == 0)
  {
    /* virgo */
    snprintf(channel, LAL_AGGREGATION_CHANNEL_MAX, "%s:%s", ifo, \
      LAL_ONLINE_VIRGO_STRAIN_CHANNEL);
  }
  else if ((strncmp(ifo, "H1", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "H2", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "L1", LAL_AGGREGATION_IFO_MAX) == 0))
  {
    /* ligo */
    snprintf(channel, LAL_AGGREGATION_CHANNEL_MAX, "%s:%s", ifo, \
      LAL_ONLINE_STRAIN_CHANNEL);
  }
  else
  {
    /* unsupported ifo */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get strain data time series */
  series = XLALFrStreamReadREAL8TimeSeries(stream, channel, start, duration, 0);
  if (series == NULL)
  {
    /* failed to read data */
    XLALFrStreamClose(stream);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* close stream */
  XLALFrStreamClose(stream);

  return series;
}


/* return data quality vector time series for given ifo, gps time and
 * duration */
INT4TimeSeries *XLALAggregationDQVector(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LALFrStream *stream;
  INT4TimeSeries *series;
  CHAR channel[LAL_AGGREGATION_CHANNEL_MAX];
  LIGOTimeGPS time_now;
  LIGOTimeGPS *latest;
  INT4 end_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get current gps time */
  if (XLALGPSTimeNow(&time_now) == NULL)
  {
    /* failed to get current time */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* check that requested data is not in the future */
  if (XLALGPSCmp(&time_now, start) == -1)
  {
    /* requested time in the future */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* determine gps time of latest frame file written */
  latest = XLALAggregationLatestGPS(ifo);
  if (latest == NULL)
  {
    /* failed to determine gps time */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* get end time of requested data */
  end_time = start->gpsSeconds + (INT4)round(duration);

  /* check that requested data has been written */
  if (latest->gpsSeconds < end_time)
  {
    /* requested data has not been written yet */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

    /* get channel name */
  if (strncmp(ifo, "V1", LAL_AGGREGATION_IFO_MAX) == 0)
  {
    /* virgo */
    snprintf(channel, LAL_AGGREGATION_CHANNEL_MAX, "%s:%s", ifo, \
      LAL_ONLINE_VIRGO_DQ_VECTOR);
  }
  else if ((strncmp(ifo, "H1", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "H2", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "L1", LAL_AGGREGATION_IFO_MAX) == 0))
  {
    /* ligo */
    snprintf(channel, LAL_AGGREGATION_CHANNEL_MAX, "%s:%s", ifo, \
      LAL_ONLINE_DQ_VECTOR);
  }
  else
  {
    /* unsupported ifo */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get data quality vector time series */
  series = XLALFrStreamReadINT4TimeSeries(stream, channel, start, duration, 0);
  if (series == NULL)
  {
    /* failed to read data */
    XLALFrStreamClose(stream);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* close stream */
  XLALFrStreamClose(stream);

  return series;
}


/* return state vector time series for given ifo, gps time, and duration */
INT4TimeSeries *XLALAggregationStateVector(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LALFrStream *stream;
  REAL4TimeSeries *state;
  INT4TimeSeries *series;
  CHAR channel[LAL_AGGREGATION_CHANNEL_MAX];
  UINT4 i;
  LIGOTimeGPS time_now;
  LIGOTimeGPS *latest;
  INT4 end_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get current gps time */
  if (XLALGPSTimeNow(&time_now) == NULL)
  {
    /* failed to get current time */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* check that requested data is not in the future */
  if (XLALGPSCmp(&time_now, start) == -1)
  {
    /* requested time in the future */
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* determine gps time of latest frame file written */
  latest = XLALAggregationLatestGPS(ifo);
  if (latest == NULL)
  {
    /* failed to determine gps time */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* get end time of requested data */
  end_time = start->gpsSeconds + (INT4)round(duration);

  /* check that requested data has been written */
  if (latest->gpsSeconds < end_time)
  {
    /* requested data has not been written yet */
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
  {
    /* failed to open stream */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get channel name */
  if ((strncmp(ifo, "H1", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "H2", LAL_AGGREGATION_IFO_MAX) == 0) || \
      (strncmp(ifo, "L1", LAL_AGGREGATION_IFO_MAX) == 0))
  {
    /* ligo */
    snprintf(channel, LAL_AGGREGATION_CHANNEL_MAX, "%s:%s", ifo, \
      LAL_ONLINE_STATE_VECTOR);
  }
  else
  {
    /* unsupported ifo */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get state vector time series */
  state = XLALFrStreamReadREAL4TimeSeries(stream, channel, start, duration, 0);
  if (state == NULL)
  {
    /* failed to read data */
    XLALFrStreamClose(stream);
    XLAL_ERROR_NULL(XLAL_EINVAL);
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
  XLALFrStreamClose(stream);

  return series;
}


/* return strain data time series for given ifo, gps time, duration, and
 * data quality vector bitmask */
REAL8TimeSeries *XLALAggregationDQStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    INT4 dq_bitmask)
{
  /* declare variables */
  INT4TimeSeries *dq_vector;
  REAL8TimeSeries *series;
  UINT4 i;

  /* checkout arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* get data quality vector */
  dq_vector = XLALAggregationDQVector(ifo, start, duration);
  if (dq_vector == NULL)
  {
    /* failed to get data quality vector */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* check for required bitmask */
  for (i = 0; i < dq_vector->data->length; i++)
  {
    if ((dq_vector->data->data[i] & dq_bitmask) != dq_bitmask)
    {
      /* invalid bitmask */
      XLALDestroyINT4TimeSeries(dq_vector);
      XLAL_ERROR_NULL(XLAL_EINVAL);
    }
  }

  /* destroy data quality vector */
  XLALDestroyINT4TimeSeries(dq_vector);

  /* get strain data time series */
  series = XLALAggregationStrainData(ifo, start, duration);
  if (series == NULL)
  {
    /* failed to get strain data time series */
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  return series;
}


/* return start position of data gap */
UINT4 XLALAggregationDQGapStart(INT4TimeSeries *series,
    INT4 dq_bitmask)
{
  /* declare variables */
  UINT4 i;
  UINT4 gap = 0;

  /* check arguments */
  if (!series)
    XLAL_ERROR(XLAL_EFAULT);

  /* check for required bitmask */
  for (i = 0; i < series->data->length; i++)
  {
    if ((series->data->data[i] & dq_bitmask) == dq_bitmask)
    {
      /* data matches bitmask */
      gap = i;
    }
    else
    {
      /* bad data */
      continue;
    }
  }

  /* return gap */
  if (gap != 0)
  {
    gap += 1;
    return gap;
  }
  else
    return 0;
}


/* return end position of data gap */
UINT4 XLALAggregationDQGapEnd(INT4TimeSeries *series,
    INT4 dq_bitmask)
{
  /* declare variables */
  UINT4 i;
  UINT4 gap = 0;

  /* check arguments */
  if (!series)
    XLAL_ERROR(XLAL_EFAULT);

  /* check for required bitmask */
  for (i = 0; i < series->data->length; i++)
  {
    if ((series->data->data[i] & dq_bitmask) == dq_bitmask)
    {
      /* data matches bitmask */
      continue;
    }
    else
    {
      /* bad data */
      gap = i;
    }
  }

  /* return gap */
  if (gap != 0)
  {
    gap += 1;
    return gap;
  }
  else
    return 0;
}


/* return end position of data gap - deprecated */
UINT4 XLALAggregationDQGap(INT4TimeSeries *series,
    INT4 dq_bitmask)
{
  /* declare variables */
  UINT4 gap = 0;

  /* check arguments */
  if (!series)
    XLAL_ERROR(XLAL_EFAULT);

  /* deprecation warning */
  XLALPrintDeprecationWarning("XLALAggregationDQGap", "XLALAggregationDQGapEnd");

  /* get end of data gap */
  gap = XLALAggregationDQGapEnd(series, dq_bitmask);

  return gap;
}


/* return strain data time series for given ifo, gps time, duration, and
 * a maximum wait time */
REAL8TimeSeries *XLALAggregationStrainDataWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait)
{
  /* declare variables */
  LALFrStream *stream;
  REAL8TimeSeries *series;
  UINT4 wait_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
    XLAL_ERROR_NULL(XLAL_EIO);

  /* initialise wait_time */
  wait_time = 0;
  do
  {
    /* try to read data */
    series = XLALAggregationStrainData(ifo, start, duration);
    if ((series == NULL) && (wait_time > max_wait))
    {
      /* already waited for maximum duration */
      XLALFrStreamClose(stream);
      XLAL_ERROR_NULL(XLAL_EIO);
    }
    else if (series == NULL)
    {
      /* failed to get series, wait */
      wait_time += LAL_ONLINE_FRAME_DURATION;
      sleep(LAL_ONLINE_FRAME_DURATION);
    }
  } while (series == NULL);

  /* close frame stream */
  XLALFrStreamClose(stream);

  return series;
}


/* return data quality vector time series for given ifo, gps time,
 * duration, and a maximum wait time */
INT4TimeSeries *XLALAggregationDQVectorWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait)
{
  /* declare variables */
  LALFrStream *stream;
  INT4TimeSeries *series;
  UINT4 wait_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
    XLAL_ERROR_NULL(XLAL_EIO);

  /* initialise wait_time */
  wait_time = 0;
  do
  {
    /* try to read data */
    series = XLALAggregationDQVector(ifo, start, duration);
    if ((series == NULL) && (wait_time > max_wait))
    {
      /* already waited for maximum duration */
      XLALFrStreamClose(stream);
      XLAL_ERROR_NULL(XLAL_EIO);
    }
    else if (series == NULL)
    {
      /* failed to get series, wait */
      wait_time += LAL_ONLINE_FRAME_DURATION;
      sleep(LAL_ONLINE_FRAME_DURATION);
    }
  } while (series == NULL);

  /* close frame stream */
  XLALFrStreamClose(stream);

  return series;
}


/* return state vector time series for given ifo, gps time, duration,
 * and a maximum wait time */
INT4TimeSeries *XLALAggregationStateVectorWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait)
{
  /* declare variables */
  LALFrStream *stream;
  INT4TimeSeries *series;
  UINT4 wait_time;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR_NULL(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR_NULL(XLAL_EFAULT);

  /* open frame stream */
  stream = XLALAggregationFrameStream(ifo, start, duration);
  if (stream == NULL)
    XLAL_ERROR_NULL(XLAL_EIO);

  /* initialise wait_time */
  wait_time = 0;
  do
  {
    /* try to read data */
    series = XLALAggregationStateVector(ifo, start, duration);
    if ((series == NULL) && (wait_time > max_wait))
    {
      /* already waited for maximum duration */
      XLALFrStreamClose(stream);
      XLAL_ERROR_NULL(XLAL_EIO);
    }
    else if (series == NULL)
    {
      /* failed to get series, wait */
      wait_time += LAL_ONLINE_FRAME_DURATION;
      sleep(LAL_ONLINE_FRAME_DURATION);
    }
  } while (series == NULL);

  /* close frame stream */
  XLALFrStreamClose(stream);

  return series;
}


/* check that all frames files, for requested data segment, are
 * available */
INT4 XLALAggregationStatFiles(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration)
{
  /* declare variables */
  LIGOTimeGPS time_now;
  LALCache *cache;
  UINT4 i;

  /* check arguments */
  if (!ifo)
    XLAL_ERROR(XLAL_EFAULT);
  if (!start)
    XLAL_ERROR(XLAL_EFAULT);

  /* get current gps time */
  if (XLALGPSTimeNow(&time_now) == NULL)
  {
    /* failed to get current time */
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* check that requested data is not in the future */
  if (XLALGPSCmp(&time_now, start) == -1)
  {
    /* requested time in the future */
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* generate frame cache for requested data */
  cache = XLALAggregationFrameCache(ifo, start, duration);
  if (cache == NULL)
  {
    /* failed to get cache */
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* loop through files in cache */
  for (i = 0; i < cache->length; i++)
  {
    /* declare variables */
    struct stat file_status;
    CHAR *filename;

    /* strip file://localhost from url */
    filename = cache->list[i].url + 16;

    /* check that file exists */
    if (stat(filename, &file_status) == -1)
    {
      /* file doesn't exist */
      XLAL_ERROR(XLAL_EIO);
    }
  }

  /* close cache */
  XLALDestroyCache(cache);

  return 0;
}
