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

/* return frame start for given gps time */
static INT4 return_frame_start(LIGOTimeGPS *gps)
{
  /* declare variables */
  static INT4 frame_start;

  /* frame start must be a multiple of ONLINE_FRAME_DURATION */
  frame_start = (INT4)floor(gps->gpsSeconds / ONLINE_FRAME_DURATION) * \
                ONLINE_FRAME_DURATION;

  return frame_start;
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


/* return path for given GPS time */
CHAR *XLALAggregationDirectoryPath(CHAR *ifo, LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationDirectoryPath";

  /* declare variables */
  CHAR *base_dir;
  CHAR *type;
  INT4 gps_dir, frame_start;
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

  /* determine gps directory name, frame start must be multiple of 16 */
  frame_start = return_frame_start(gps);
  gps_dir = (INT4)floor(frame_start / 100000);

  /* construct directory */
  LALSnprintf(directory, FILENAME_MAX, "%s/%c-%s-%d", base_dir, ifo[0], \
      type, gps_dir);

  return directory;
}


/* return frame filename for given gps time and ifo */
CHAR *XLALAggregationFrameFilename(CHAR *ifo, LIGOTimeGPS *gps)
{
  static const char *func = "XLALAggregationFrameFilename";

  /* declare variables */
  CHAR *type;
  INT4 frame_start;
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
  frame_start = return_frame_start(gps);

  /* construct frame filename */
  LALSnprintf(filename, FILENAME_MAX, "%c-%s-%d-%d.gwf", ifo[0], type, \
      frame_start, ONLINE_FRAME_DURATION);

  return filename;
}


/* return full path to frame for a given gps time and ifo */
CHAR *XLALAggregationFramePathFilename(CHAR *ifo, LIGOTimeGPS *gps)
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


/* return url to frame for a given gps time and ifo */
CHAR *XLALAggregationFrameURL(CHAR *ifo, LIGOTimeGPS *gps)
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


/* return frame cache of required frames */
FrCache *XLALAggregationFrameCache(CHAR *ifo, LIGOTimeGPS *start, INT4 length)
{
  static const char *func = "XLALAggregationFrameCache";

  /* declare variables */
  LIGOTimeGPS gps;
  CHAR *type;
  CHAR *url;
  FrCache *cache;
  int i = 0;
  INT4 frame_start;
  INT4 last_frame_start;
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
  frame_start = return_frame_start(start);
  last_frame_start = return_frame_start(&gps);
  frame_duration = (last_frame_start + ONLINE_FRAME_DURATION) - frame_start;
  num_frames = frame_duration / ONLINE_FRAME_DURATION;

  /* initilise cache */
  cache = LALCalloc(1, sizeof(*cache));
  cache->numFrameFiles = num_frames;
  cache->frameFiles = LALCalloc(num_frames, sizeof(*cache->frameFiles));
  if (!cache->frameFiles)
  {
    LALFree(cache);
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
    frame_start = return_frame_start(&gps);

    /* allocate memory for cache entry */
    file->source = LALMalloc(strlen(ifo) + 1);
    if (!file->source)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }
    file->description = LALMalloc(strlen(type) + 1);
    if (!file->description)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }
    file->url = LALMalloc(strlen(url) + 1);
    if (!file->url)
    {
      XLALFrDestroyCache(cache);
      XLAL_ERROR_NULL(func, XLAL_ENOMEM);
    }

    /* add frame to cache */
    strcpy(file->source, ifo);
    strcpy(file->description, type);
    file->startTime = frame_start;
    file->duration = ONLINE_FRAME_DURATION;
    strcpy(file->url, url);
  }

  return cache;
}


/* return required frame stream */
FrStream *XLALAggregationFrameStream(CHAR *ifo, LIGOTimeGPS *start, INT4 length)
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
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* destroy cache */
  XLALFrDestroyCache(cache);

  return stream;
}


/* return strain data time series for given ifo, start time, and duration */
REAL8TimeSeries *XLALAggregationStrainData(CHAR *ifo, LIGOTimeGPS *start, INT4 length)
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
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "%s:%s", ifo, ONLINE_STRAIN_CHANNEL);
  series = XLALFrReadREAL8TimeSeries(stream, channel, start, (REAL8)length, num_points);
  if (series == NULL)
  {
    /* failed to read data */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* close stream */
  XLALFrClose(stream);

  return series;
}


/* return data quality vector time series for given ifo, start time, and duration */
INT4TimeSeries *XLALAggregationDQVector(CHAR *ifo, LIGOTimeGPS *start, INT4 length)
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

  /* get strain data time series */
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "%s:%s", ifo, ONLINE_DQ_VECTOR);
  series = XLALFrReadINT4TimeSeries(stream, channel, start, (REAL8)length, num_points);
  if (series == NULL)
  {
    /* failed to read data */
    XLAL_ERROR_NULL(func, XLAL_EINVAL);
  }

  /* close stream */
  XLALFrClose(stream);

  return series;
}
