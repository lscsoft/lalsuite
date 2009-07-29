/*
 * online_datafind.c - online frame cache generator
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
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/Aggregation.h>
#include <lal/XLALError.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/Date.h>

#include <lalapps.h>

/* flags for getopt_long */
extern int vrbflg;

/* global variables */
CHAR *ifo = NULL;
LIGOTimeGPS gps_start = {0, 0};
LIGOTimeGPS gps_end = {0, 0};
INT4 timeout = 0;
CHAR *output_filename = NULL;
CHAR *observatory = NULL;
CHAR *frame_type = NULL;

/*
 * helper functions
 */

/* parse command line options */
static void parse_options(INT4 argc, CHAR *argv[])
{
  int c = -1;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"verbose", no_argument, &vrbflg, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'a'},
      {"debug-level", required_argument, 0, 'b'},
      {"ifo", required_argument, 0, 'c'},
      {"gps-start-time", required_argument, 0, 'd'},
      {"gps-end-time", required_argument, 0, 'e'},
      {"timeout", required_argument, 0, 'f'},
      {"output", required_argument, 0, 'g'},
      {"observatory", required_argument, 0, 'h'},
      {"type", required_argument, 0, 'i'},
      {"url-type", required_argument, 0, 'j'},
      {"server", required_argument, 0, 'k'},
      {"lal-cache", no_argument, 0, 'l'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    /* parse options */
    c = getopt_long_only(argc, argv, "ab:c:d:e:f:g:h:i:j:k:l", \
        long_options, &option_index);

    if (c == -1)
    {
      /* end of options, break loop */
      break;
    }

    switch(c)
    {
      case 0:
        /* if this option sets a flag, do nothing else now */
        if (long_options[option_index].flag != 0)
        {
          break;
        }
        else
        {
          fprintf(stderr, "error parsing option %s with argument %s\n", \
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'a':
        /* help */
        fprintf(stdout, "Usage: lalapps_online_datafind [options]\n");
        fprintf(stdout, " --help                 print this message\n");
        fprintf(stdout, " --verbose              run in verbose mode\n");
        fprintf(stdout, " --debug-level N        set lalDebugLevel\n");
        fprintf(stdout, " --ifo IFO              set IFO\n");
        fprintf(stdout, " --gps-start-time GPS   start of GPS time range\n");
        fprintf(stdout, " --gps-end-time GPS     end of GPS time range\n");
        fprintf(stdout, " --timeout TIME         set timeout\n");
        fprintf(stdout, " --output FILE          output to FILE\n");
        fprintf(stdout, " --observatory SITE     set SITE\n");
        fprintf(stdout, " --type FRAME_TYPE      set FRAME_TYPE\n");
        exit(0);
        break;

      case 'b':
        /* get debug level */
        lalDebugLevel = atoi(optarg);
        break;

      case 'c':
        /* get ifo */
        optarg_len = strlen(optarg) + 1;
        ifo = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifo, optarg, optarg_len);
        break;

      case 'd':
        /* get gps start time */
        gps_start.gpsSeconds = atoi(optarg);
        gps_start.gpsNanoSeconds = 0;
        if (gps_start.gpsSeconds <= 0)
        {
          fprintf(stderr, "invalid argument to --%s: %d\n", \
              long_options[option_index].name, gps_start.gpsSeconds);
          exit(1);
        }
        break;

      case 'e':
        /* get gps end time */
        gps_end.gpsSeconds = atoi(optarg);
        gps_end.gpsNanoSeconds = 0;
        if (gps_end.gpsSeconds <= 0)
        {
          fprintf(stderr, "invalid argument to --%s: %d\n", \
              long_options[option_index].name, gps_end.gpsSeconds);
          exit(1);
        }
        break;

      case 'f':
        /* get timeout */
        timeout = atoi(optarg);
        if (timeout < 0)
        {
          fprintf(stderr, "invalid argument to --%s: %d\n", \
              long_options[option_index].name, timeout);
          exit(1);
        }

      case 'g':
        /* get filename */
        optarg_len = strlen(optarg) + 1;
        output_filename = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(output_filename, optarg, optarg_len);
        break;

      case 'h':
        /* set observatory */
        optarg_len = strlen(optarg) + 1;
        observatory = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(observatory, optarg, optarg_len);
        break;

      case 'i':
        /* set frame type */
        optarg_len = strlen(optarg) + 1;
        frame_type = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(frame_type, optarg, optarg_len);
        break;

      case 'j':
      case 'k':
      case 'l':
        /* url type, server, and lal-cache have no effect, they are just
         * provided for compatibility with ligo_data_find */
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "unknown error while parsing options\n");
        exit(1);
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "extraneous command line arguments:\n");
    while(optind < argc)
    {
      fprintf(stderr, "%s\n", argv[optind++]);
    }
    exit(1);
  }

  /*
   * check for required arguments
   */

  /* are observatory and type both specifed */
  if (((observatory != NULL) && (frame_type == NULL)) || \
      ((observatory == NULL) && (frame_type != NULL)))
  {
    fprintf(stdout, "both --observatory and --type must be specified\n");
    exit(1);
  }

  /* are ifo and observatory both specified */
  if ((ifo != NULL) && (observatory != NULL))
  {
    fprintf(stdout, "--ifo and --observatory are both specified\n");
    exit(1);
  }

  /* ifo or observatory */
  if ((ifo == NULL) && (observatory == NULL))
  {
    fprintf(stderr, "--ifo or --observatory and --type must be specified\n");
    exit(1);
  }

  /* gps start time */
  if (gps_start.gpsSeconds == 0)
  {
    fprintf(stderr, "--gps-start-time must be specified\n");
    exit(1);
  }

  /* gps end time */
  if (gps_end.gpsSeconds == 0)
  {
    fprintf(stderr, "--gps-end-time must be specified\n");
    exit(1);
  }

  return;
}


/* main program entry point */
INT4 main(INT4 argc, CHAR *argv[])
{
  /* declare variables */
  LIGOTimeGPS *latest_time;
  LIGOTimeGPS time_now;
  INT4 duration;
  FrCache *cache;
  CHAR *type;
  INT4 wait_time;
  CHAR *ptimeout;
  INT4 total_wait = 0;

  /* get maximum wait time from ONLINEHOFT_TIMEOUT */
  ptimeout = getenv("ONLINEHOFT_TIMEOUT");
  if (ptimeout != NULL)
  {
    /* get timout from environment */
    timeout = atoi(ptimeout);
  }

  /* parse command line options */
  parse_options(argc, argv);

  /* get duration */
  duration = (INT4)round(XLALGPSDiff(&gps_end, &gps_start));

  /* get ifo, if observatory is specified */
  if (ifo == NULL)
  {
    /* allocate memory for ifo */
    ifo = (CHAR *)calloc(LIGOMETA_IFO_MAX, sizeof(CHAR));

    /* get ifo from frame_type */
    strncpy(ifo, frame_type, 2);
  }

  /* get time of gps time of latest frame */
  latest_time = XLALAggregationLatestGPS(ifo);
  if (latest_time == NULL)
  {
    fprintf(stderr, "error: unable to determine latest GPS time\n");
    exit(1);
  }

  /* get current gps time */
  if (XLALGPSTimeNow(&time_now) == NULL)
  {
    fprintf(stderr, "error: unable to determine current time\n");
    exit(1);
  }

  /* report time info */
  if (vrbflg)
  {
    fprintf(stdout, "current time:          %d\n", time_now.gpsSeconds);
    fprintf(stdout, "latest data available: %d\n", latest_time->gpsSeconds);
    fprintf(stdout, "requested start:       %d\n", gps_start.gpsSeconds);
    fprintf(stdout, "requested end:         %d\n", gps_end.gpsSeconds);
    fprintf(stdout, "requested duration:    %9d\n", duration);
  }

  /* is requested data in the future? */
  if (XLALGPSCmp(&time_now, &gps_end) == -1)
  {
    /* determine wait time */
    wait_time = (INT4)round(XLALGPSDiff(&gps_end, &time_now));

    /* wait for data to be available */
    fprintf(stdout, "requested data is in the future, waiting: %ds\n", wait_time);
    sleep(wait_time);
  }

  /* has requested data been written to disk */
  if (XLALGPSCmp(latest_time, &gps_end) == -1)
  {
    /* determine wait time */
    wait_time = (INT4)round(XLALGPSDiff(&gps_end, latest_time));

    /* does wait exceed timeout? */
    if (wait_time > timeout)
    {
      fprintf(stderr, "data unavailable, wait exceeds timeout: %ds\n", wait_time);
      exit(1);
    }

    /* wait for data to be available */
    fprintf(stdout, "data unavailable, waiting %ds\n", wait_time);
    sleep(wait_time);
    total_wait += wait_time;

    /* loop to wait for requested data */
    do
    {
      /* get new latest gps time */
      LIGOTimeGPS *latest_gps;
      latest_gps = XLALAggregationLatestGPS(ifo);
      if (latest_gps == NULL)
      {
        fprintf(stderr, "error: unable to determine latest GPS time\n");
        exit(1);
      }

      /* has requested data been written to disk? */
      if (XLALGPSCmp(latest_gps, &gps_end) == -1)
      {
        /* determine wait time */
        wait_time = (INT4)round(XLALGPSDiff(&gps_end, latest_gps));

        /* does required wait exceed timeout? */
        if ((total_wait + wait_time) > timeout)
        {
          fprintf(stderr, "data unavailable, wait exceeds timeout: %ds\n", \
              total_wait + wait_time);
          exit(1);
        }

        /* wait for data to be available */
        fprintf(stdout, "data unavailable, waiting %ds\n", wait_time);
        sleep(wait_time);
        total_wait += wait_time;
      }
      else
      {
        /* check that all files have been transferred, stat frames */
        if (XLALAggregationStatFiles(ifo, &gps_start, duration) == 0)
        {
          /* data is available, break do-while loop */
          break;
        }
      }
    } while(total_wait < timeout);
  }

  /* get frame cache */
  cache = XLALAggregationFrameCache(ifo, &gps_start, duration);
  if (cache == NULL)
  {
    fprintf(stderr, "error: failed to get frame cache\n");

    /* free memory and exit */
    if (ifo)
      free(ifo);
    if (observatory)
      free(observatory);
    if (frame_type)
      free(frame_type);
    exit(xlalErrno);
  }

  /* get frame type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    fprintf(stderr, "error: failed to get frame type\n");

    /* free memory and exit */
    if (ifo)
      free(ifo);
    if (observatory)
      free(observatory);
    if (frame_type)
      free(frame_type);
    exit(xlalErrno);
  }

  if (output_filename == NULL)
  {
    /* allocate memory for filename */
    output_filename = (CHAR *)calloc(FILENAME_MAX, sizeof(CHAR));

    /* create name for cache file */
    snprintf(output_filename, FILENAME_MAX, "%c-%s-%d-%d.cache", \
        ifo[0], type, gps_start.gpsSeconds, gps_end.gpsSeconds);
    fprintf(stdout, "output filename not specified, using: %s\n", \
        output_filename);
  }

  /* save cache */
  XLALFrExportCache(cache, output_filename);

  /* free memory */
  if (ifo)
    free(ifo);
  if (observatory)
    free(observatory);
  if (frame_type)
    free(frame_type);
  free(output_filename);

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  return 0;
}
