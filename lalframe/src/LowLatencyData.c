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
 * Copyright (C) 2011 Adam Mercer, Leo Singer
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <fnmatch.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#if HAVE_INOTIFY
#include <sys/inotify.h>
#else
#include <errno.h>
#include <dirent.h>
#endif

#include <lal/LowLatencyData.h>

/**
   \ingroup LowLatencyData_h
*/

struct tagLowLatencyData
{
  /* File descriptor for inotify or directory */
  int fd;

#if HAVE_INOTIFY
  /* Watch descriptor for inotify */
  int wd;
  /* Buffer for change notifications */
  char buf[sizeof(struct inotify_event) + FILENAME_MAX + 1];
  /* Pointer to next unconsumed inotify event in buf */
  struct inotify_event *next;
#else
  DIR *dir;
#endif

  /* Path to directory */
  char path[FILENAME_MAX + 1];
  /* Filename pattern */
  char pattern[FILENAME_MAX + 1];
  /* Last filename successfully read */
  char last_filename[FILENAME_MAX + 1];

  /* Pointer frame in memory, if any */
  void *frame_ptr;
  /* Size of frame in memory */
  size_t frame_size;
};


/**
 * Create a new instance of LowLatencyData.
 * It will be ready to read the next available frame file.
 */
LowLatencyData *XLALLowLatencyDataOpen(
  const char *data_path,    /* for example:  "/dev/shm/H1"    */
  const char *observatory,  /* for example:  "H"              */
  const char *frame_type    /* for example:  "H1_DMT_C00_L0"  */
)
{
  LowLatencyData *data = (LowLatencyData *)calloc(1, sizeof(LowLatencyData));
  if (!data)
    return NULL;

  strcpy(data->path, data_path);

  snprintf(data->pattern, FILENAME_MAX, "%s-%s-*-*.gwf", observatory, frame_type);

#if HAVE_INOTIFY
  data->fd = inotify_init();
  if (data->fd == -1)
  {
    free(data);
    return NULL;
  }

  data->wd = inotify_add_watch(data->fd, data_path, IN_MOVED_TO);
  if (data->wd == -1)
  {
    close(data->fd);
    free(data);
    return NULL;
  }

  data->next = NULL;
#else
  data->dir = opendir(data_path);
  if (!data->dir)
  {
    free(data);
    return NULL;
  }

  data->fd = dirfd(data->dir);
  if (data->fd == -1)
  {
    closedir(data->dir);
    free(data);
    return NULL;
  }
#endif

  data->frame_ptr = MAP_FAILED;

  /* technically, doesn't have to be set because calloc zeros the struct,
   * but best to make this explicit */
  data->last_filename[0] = '\0';

  return data;
}


/**
 * Destroy an instance of LowLatencyData.
 * Release resources associated with monitoring /dev/shm.
 */
void XLALLowLatencyDataClose(LowLatencyData *data)
{
#if HAVE_INOTIFY
  inotify_rm_watch(data->fd, data->wd);
  data->wd = -1;
  close(data->fd);
  data->fd = -1;
#else
  closedir(data->dir);
  data->dir = NULL;
  data->fd = -1;
#endif
  if (data->frame_ptr != MAP_FAILED)
  {
    munmap(data->frame_ptr, data->frame_size);
    data->frame_ptr = MAP_FAILED;
  }
  free(data);
}


static char *XLALLowLatencyDataNextFilename(LowLatencyData *data)
{
#if HAVE_INOTIFY /* Inotify version */

  char *filename;

  do {
    int match;
    struct inotify_event *evt;

    /* Ignore queue overflow events. */
    do {

      /* Check to see if there is another inotify event already in our buffer.
       * If so, return that event and advance to the next event.
       * Otherwise, try to read from the inotify file descriptor again.*/

      if (data->next && (char *) data->next < (char *) &data->buf + sizeof(data->buf))
        evt = data->next;
      else if (read(data->fd, data->buf, sizeof(data->buf)) != -1)
        evt = (struct inotify_event *) data->buf;
      else
        return NULL; /* Something went wrong reading from inotify. */

      data->next = (struct inotify_event *) ((char *) evt + sizeof(struct inotify_event) + evt->len);

    } while (!(evt->mask & ~IN_Q_OVERFLOW));

    match = fnmatch(data->pattern, evt->name, 0);
    if (match == 0) /* match found */
      filename = evt->name;
    else if (match == FNM_NOMATCH) /* not a match */
      filename = NULL;
    else /* error occurred in fnmatch */
      return NULL;

  } while (!filename || (data->last_filename && strcmp(filename, data->last_filename) <= 0));

#else /* Directory polling version */

  char filename[FILENAME_MAX + 1];

  do /* Loop until candidate filename is a non-empty string. */ {
    struct dirent *entry;

    for (filename[0] = '\0', rewinddir(data->dir); errno = 0, entry = readdir(data->dir); )
    {
      int match;

      /* Demand that the filename matches the pattern. */
      match = fnmatch(data->pattern, entry->d_name, 0);
      if (match == FNM_NOMATCH) /* not a match */
        continue;
      else if (match != 0) /* error occurred in fnmatch */
        return NULL;

      /* Demand that the filename is lexicographically greater than the
       * last opened filename. */
      if (strcmp(entry->d_name, data->last_filename) < 1)
        continue;

      /* If there is already another candidate filename, then we have to
       * compare the name of this entry lexicographically with it. */
      if (filename[0])
      {
        if (data->last_filename[0])
        {
          /* If we have returned a filename before, then demand that
           * the new filename is lexicographically less than the candidate. */
          if (strcmp(entry->d_name, filename) >= 0)
            continue;
        } else {
          /* Otherwise, we are searching for the first filename of the
           * session, so we want the newest possible file, so we demand
           * that the new filename is lexographically greater than the
           * candidate. */
          if (strcmp(entry->d_name, filename) <= 0)
            continue;
        }
      }

      /* This entry is even better than the previous candidate, so make
       * it the new candidate. */
      strcpy(filename, entry->d_name);
    }

    /* If an error occurred when calling readdir, fail. */
    if (errno)
      return NULL;

    /* If no new filename was found, then sleep 1 second before returning. */
    if (!filename[0])
      sleep(1);

  } while (!filename[0]);

#endif

  /* Record this as the last filename. */
  strcpy(data->last_filename, filename);

  /* Success! */
  return data->last_filename;
}


/**
 * Retrieve the next available frame file.
 * On success, a pointer to a new in-memory frame file is returned.  If size is
 * a non-NULL pointer, then the size of the file in bytes is written to (*size).
 * Any prevoiusly returned pointer from this instance of LowLatencyData is
 * invalidated.
 *
 * On failure, NULL is returned and (*size) is not modified.  It is unspecified
 * whether any previously returned pointer is invalidated on failure.
 */
void *XLALLowLatencyDataNextBuffer(LowLatencyData *data, size_t *size)
{
  int fd;
  struct stat st;
  char pathname[FILENAME_MAX + 1];
  char *filename = XLALLowLatencyDataNextFilename(data);

  if (!filename)
    return NULL;

  /* create full path to file */
  {
    int ret = snprintf(pathname, sizeof(pathname), "%s/%s", data->path, filename);
    if (ret < 0 || (size_t) ret >= (size_t) sizeof(pathname))
      return NULL;
  }

  /* open file */
  fd = open(pathname, O_RDONLY);
  if (fd == -1) return NULL;

  /* get stats */
  if (fstat(fd, &st) == -1)
  {
    close(fd);
    return NULL;
  }

  /* unmap old file */
  if (data->frame_ptr != MAP_FAILED)
  {
    munmap(data->frame_ptr, data->frame_size);
    data->frame_ptr = MAP_FAILED;
  }

  data->frame_size = st.st_size;
  data->frame_ptr = mmap(0, data->frame_size, PROT_READ, MAP_SHARED, fd, 0);
  close(fd);

  if (data->frame_ptr == MAP_FAILED)
    return NULL;

  if (size)
    *size = data->frame_size;

  return data->frame_ptr;
}
