/*
 * lld.c - low latency data routine prototype code
 *       - functions
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
 * Copyright (C) 2011 Adam Mercer
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>
#include <getopt.h>
#include <math.h>

#include <lal/LowLatencyData.h>

/* default options */
#define LLD_DATA_PATH "/dev/shm"
#define LLD_OBSERVATORY "X"
#define LLD_FRAME_TYPE "R"

/* global variables */
int vrbflg;
LIGOTimeGPS requested_time = {0,0};
CHAR *observatory = NULL;
CHAR *frame_type = NULL;
CHAR *data_path = NULL;

/*
 * local helper functions
 */

/* prototypes */
void parse_options(int argc, char **argv);
void display_help_message(void);

/* parse command line options */
void parse_options(
    int argc,
    char **argv)
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
      {"time", required_argument, 0, 'b'},
      {"observatory", required_argument, 0, 'c'},
      {"type", required_argument, 0, 'd'},
      {"path", required_argument, 0, 'e'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    /* parse options */
    c = getopt_long_only(argc, argv, \
        "ab:c:d:", \
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
          break;
        else
        {
          fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'a':
        /* help! */
        display_help_message();
        exit(0);
        break;

      case 'b':
        /* time */
        requested_time.gpsSeconds = atoi(optarg);
        if (requested_time.gpsSeconds <= 0)
        {
          fprintf(stderr, "invalid argument to --%s: %d\n", \
              long_options[option_index].name, requested_time.gpsSeconds);
          exit(1);
        }
        break;

      case 'c':
        /* observatory */
        optarg_len = strlen(optarg) + 1;
        observatory = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(observatory, optarg, optarg_len);
        break;

      case 'd':
        /* frame type */
        optarg_len = strlen(optarg) + 1;
        frame_type = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(frame_type, optarg, optarg_len);
        break;

      case 'e':
        /* data path */
        optarg_len = strlen(optarg) + 1;
        data_path = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(data_path, optarg, optarg_len);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "banana in disk drive error\n");
        exit(1);
        break;
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "extraneous command line arguments:\n");
    while(optind < argc)
      fprintf(stderr, "%s\n", argv[optind++]);
    exit(1);
  }

  /* assign defaults if not specified */

  /* observatory */
  if (observatory == NULL)
  {
    observatory = (CHAR *)calloc(strlen(LLD_OBSERVATORY) + 1, sizeof(CHAR));
    memcpy(observatory, LLD_OBSERVATORY, strlen(LLD_OBSERVATORY) + 1);
  }

  /* frame type */
  if (frame_type == NULL)
  {
    frame_type = (CHAR *)calloc(strlen(LLD_FRAME_TYPE) + 1, sizeof(CHAR));
    memcpy(frame_type, LLD_FRAME_TYPE, strlen(LLD_FRAME_TYPE) + 1);
  }

  /* data path */
  if (data_path == NULL)
  {
    data_path = (CHAR *)calloc(strlen(LLD_DATA_PATH) + 1, sizeof(CHAR));
    memcpy(data_path, LLD_DATA_PATH, strlen(LLD_DATA_PATH) + 1);
  }

  /*
   * check for required arguments
   */

  /* time */
  if (requested_time.gpsSeconds == 0)
  {
    fprintf(stderr, "--time must be specified\n");
    exit(1);
  }

  return;
}

/* display help */
void display_help_message(void)
{
  fprintf(stdout, "Usage: lalframe_lld [options]\n");
  fprintf(stdout, " --help          print this message\n");
  fprintf(stdout, " --verbose       run in verbose mode\n");
  fprintf(stdout, " --time          time to return\n");
  fprintf(stdout, " --observatory   observatory [default = X]\n");
  fprintf(stdout, " --type          frame type [default = R]\n");
  fprintf(stdout, " --path          location of data [default = /dev/shm]\n");
  return;
}

/*
 * main program entry point
 */

int main(int argc, char **argv)
{
  /* parse command line options */
  parse_options(argc, argv);

  /* get filename of next frame following requested time */
  CHAR *filename;
  filename = LALFrameLLDNextFrameName(data_path, observatory, frame_type, requested_time);
  if (filename == NULL)
  {
    fprintf(stderr, "Unable to get frame name\n");

    /* free memory */
    free(observatory);
    free(frame_type);

    /* exit */
    exit(1);
  }

  /* report */
  fprintf(stdout, "filename = %s\n", filename);

  /* free memory */
  free(observatory);
  free(frame_type);

  /* exit cleanly */
  exit(0);
}
