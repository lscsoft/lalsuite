/*
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
 * Copyright (C) 2011 Adam Mercer
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>
#include <getopt.h>
#include <math.h>

#include <lal/LALDatatypes.h>
#include <lal/LowLatencyData.h>

/*
 *  static functions
 */

/* prototypes */

/* return glob of files with given pattern */
static glob_t *lld_get_files(CHAR *glob_pattern);

/* return gps start time of first frame file */
static INT4 lld_gps_start_time(CHAR *data_path, CHAR *observatory,
    CHAR *frame_type, glob_t *files);

/* return gps start time of last frame file */
static INT4 lld_last_gps_start_time(CHAR *data_path, CHAR *observatory,
    CHAR *frame_type, glob_t *files);

/* return frame duration */
static INT4 lld_frame_duration(CHAR *data_path, CHAR *observatory,
    CHAR *frame_type, glob_t *files);

/* return gps end time of last frame file */
static INT4 lld_gps_end_time(CHAR *data_path, CHAR *observatory,
    CHAR *frame_type, glob_t *files);

/* check if requested time is in range */
static INT4 lld_time_in_range(CHAR *data_path, CHAR *observatory,
    CHAR *frame_type, INT4 requested_time, glob_t *files);

/* functions */

/* return glob of files with given pattern */
static glob_t *lld_get_files(
    CHAR *glob_pattern)
{
  /* variables */
  glob_t *glob_files;

  /* allocate memory for glob_files */
  glob_files = calloc(1, sizeof(*glob_files));

  /* get files matching pattern */
  if (glob(glob_pattern, 0, NULL, glob_files) != 0)
  {
    /* glob failed */
    globfree(glob_files);
    return NULL;
  }

  return glob_files;
}

/* return gps start time of first frame file */
static INT4 lld_gps_start_time(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    glob_t *files)
{
  /* variables */
  CHAR fmt[FILENAME_MAX];
  INT4 gps_start_time;

  /* get sscanf format for parsing start time */
  snprintf(fmt, FILENAME_MAX, "%s/%s-%s-%%d-%%*d.gwf", data_path, \
      observatory, frame_type);

  /* get start time */
  if (sscanf(files->gl_pathv[0], fmt, &gps_start_time) != 1)
  {
    /* error parsing start time */
    return -1;
  }

  return gps_start_time;
}

/* return gps start time of last frame file */
static INT4 lld_last_gps_start_time(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    glob_t *files)
{
  /* variables */
  CHAR fmt[FILENAME_MAX];
  INT4 gps_end_time;

  /* get sscanf format for parsing start time */
  snprintf(fmt, FILENAME_MAX, "%s/%s-%s-%%d-%%*d.gwf", data_path, \
      observatory, frame_type);

  /* get start time */
  if (sscanf(files->gl_pathv[files->gl_pathc - 1], fmt, &gps_end_time) != 1)
  {
    /* error parsing start time */
    return -1;
  }

  return gps_end_time;
}

/* return frame duration */
static INT4 lld_frame_duration(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    glob_t *files)
{
  /* variables */
  CHAR fmt[FILENAME_MAX];
  INT4 duration;

  /* get sscanf format for parsing start time */
  snprintf(fmt, FILENAME_MAX, "%s/%s-%s-%%*d-%%d.gwf", data_path, \
      observatory, frame_type);

  /* get frame duration */
  if (sscanf(files->gl_pathv[0], fmt, &duration) != 1)
  {
    /* error parsing frame duration */
    return -1;
  }

  return duration;
}

/* return gps end time of last frame file */
static INT4 lld_gps_end_time(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    glob_t *files)
{
  /* variables */
  INT4 last_start_time;
  INT4 duration;

  /* get start time of last frame */
  last_start_time = lld_last_gps_start_time(data_path, observatory, \
      frame_type, files);
  if (last_start_time == -1)
    return -1;

  /* get frame duration */
  duration = lld_frame_duration(data_path, observatory, frame_type, \
      files);
  if (duration == -1)
    return -1;

  return last_start_time + duration;
}

/* check if requested time is in range */
static INT4 lld_time_in_range(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    INT4 requested_time,
    glob_t *files)
{
  /* variables */
  INT4 start_time;
  INT4 end_time;

  /* get start and end times */
  start_time = lld_gps_start_time(data_path, observatory, frame_type, files);
  end_time = lld_gps_end_time(data_path, observatory, frame_type, files);
  if ((start_time == -1) || (end_time == -1))
    return -1;

  /* is requested time in range */
  if (!((start_time <= requested_time) && (end_time >= requested_time)))
    return -1;

  return 0;
}

/*
 * externally linkable functions
 */

/* return file name of the next frame following requested time */
CHAR *LALFrameLLDNextFrameName(
    CHAR *data_path,
    CHAR *observatory,
    CHAR *frame_type,
    LIGOTimeGPS *requested_time)
{
  /* variables */
  CHAR glob_pattern[FILENAME_MAX];
  glob_t *files;
  INT4 start_time;
  INT4 duration;
  REAL8 span;
  INT4 n_frames;
  INT4 frame_start;
  CHAR *filename;

  /* get glob pattern */
  snprintf(glob_pattern, FILENAME_MAX, "%s/*-*-*-*.gwf", data_path);

  /* get list of files */
  files = lld_get_files(glob_pattern);
  if (files == NULL)
    return NULL;

  /* get start time of first frame */
  start_time = lld_gps_start_time(data_path, observatory, frame_type, files);
  if (start_time == -1)
  {
    /* free memory */
    globfree(files);

    return NULL;
  }

  /* get frame duration */
  duration = lld_frame_duration(data_path, observatory, frame_type, files);
  if (duration == -1)
  {
    /* free memory */
    globfree(files);

    return NULL;
  }

  /* is requested time in range */
  if (lld_time_in_range(data_path, observatory, frame_type, requested_time->gpsSeconds, files) != 0)
  {
    /* free memory */
    globfree(files);

    return NULL;
  }

  /* determine start time of frame */
  span = requested_time->gpsSeconds - start_time;
  n_frames = (INT4)ceil(span / duration);
  frame_start = start_time + (n_frames * duration);

  /* report */
  filename = (CHAR *)calloc(FILENAME_MAX, sizeof(CHAR));
  snprintf(filename, FILENAME_MAX, "%s/%s-%s-%d-%d.gwf", data_path, observatory, frame_type, frame_start, duration);

  return filename;
}
