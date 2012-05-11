/*
 * TestLowLatencyData1.c
 * Test that the low latency data reader correctly reads a sequence of frames
 * that are written in chronological order.
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
 * Copyright (C) 2011 Leo Singer
 */


#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>

#include "config.h"
#include <lal/Date.h>
#include <lal/LowLatencyData.h>
#include <lal/LALFrameL.h>
#include <lal/LALFrameIO.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static char tempdir[] = "/tmp/TestLowLatencyDataXXXXXX";

static void write_dummy_frame(long gps_start_time, int duration)
{
  char filename[FILENAME_MAX + 1];
  LIGOTimeGPS epoch;
  FrameH *frame;

  XLALGPSSet(&epoch, gps_start_time, 0);
  sprintf(filename, "%s/X-R-%010ld-%d.gwf", tempdir, gps_start_time, duration);
  frame = XLALFrameNew(&epoch, duration, "FOOBAR", 0, 0, 0);
  XLALFrameWrite(frame, filename, 0);
}

static void *write_dummy_frames(void UNUSED *not_used)
{
  int i;
  for (i = 0; i < 1000; i ++)
    write_dummy_frame(1000000000 + i, 1);
  return NULL;
}

static double get_gps_start_time(void *bytes, size_t size)
{
  double gps;
  FrameH *frame = FrameReadFromBuf(bytes, size, 0);
  gps = frame->GTimeS;
  FrameFree(frame);
  return gps;
}

static void remove_tempdir(void)
{
  char path[FILENAME_MAX + 1];
  struct dirent *entry;
  DIR *dir = opendir(tempdir);

  if (!dir)
  {
    perror("opendir");
    return;
  }

  for (; errno = 0, entry = readdir(dir); )
  {
    if (strcmp(&entry->d_name[strlen(entry->d_name) - 4], ".gwf") == 0)
    {
      sprintf(path, "%s/%s", tempdir, entry->d_name);
      if (unlink(path) == -1)
        perror("unlink");
    }
  }

  if (errno)
    perror("readdir");

  if (closedir(dir) == -1)
    perror("closedir");

  if (rmdir(tempdir) == -1)
    perror("rmdir");
}

int main(void)
{
  pthread_t thread;
  LowLatencyData *reader;
  int i;

  /* exit() abnormally if anything goes wrong in XLAL land */
  XLALSetErrorHandler(XLALExitErrorHandler);

  if (!mkdtemp(tempdir))
  {
    perror("mkdtemp");
    exit(1);
  }

  atexit(remove_tempdir);

  reader = XLALLowLatencyDataOpen(tempdir, "X", "R");
  if (!reader)
  {
    perror("XLALLowLatencyDataOpen");
    exit(1);
  }

#if !HAVE_INOTIFY
  /* When using the readdir implementation, the first file picked up is the
   * 'newest' one currently in the directory, by contrast with the inotify
   * implementation, where the first file picked up is the first one placed
   * in the directory after the first entry into XLALLowLatencyDataNextBuffer.
   * In order for the first buffer returned by the readdir implementation to
   * be predictable, we have to write one frame and then check that the reader
   * picked it up before spawning the writer thread. */
  write_dummy_frame(999999999, 1);
  XLALLowLatencyDataNextBuffer(reader, NULL);
#endif

  pthread_create(&thread, NULL, write_dummy_frames, NULL);

  for (i = 0; i < 1000; i ++)
  {
    void *bytes;
    size_t size;
    bytes = XLALLowLatencyDataNextBuffer(reader, &size);
    if (!bytes)
    {
      perror("XLALLowLatencyDataNextBuffer");
      pthread_join(thread, NULL);
      exit(1);
    }
    double expected_gps_time = 1000000000 + i;
    double gps_time = get_gps_start_time(bytes, size);
    if (expected_gps_time != gps_time)
    {
      fprintf(stderr, "Expected GPS time %f, got GPS time %f\n", expected_gps_time, gps_time);
      pthread_join(thread, NULL);
      exit(1);
    }
  }

  XLALLowLatencyDataClose(reader);
  pthread_join(thread, NULL);

  exit(0);
}
