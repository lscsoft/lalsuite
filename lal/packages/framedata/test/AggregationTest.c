/*
 * AggregationTest.c - test online frame data aggregation routines
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

#include <lal/LALDatatypes.h>
#include <lal/Aggregation.h>
#include <lal/XLALError.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>

/* flags for getopt_long */
int vrbflg;

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
      {"help", no_argument, 0, 'a'},
      {"verbose", no_argument, &vrbflg, 1},
      {"debug-level", required_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;

    /* parse options */
    c = getopt_long_only(argc, argv, "d:", long_options, &option_index);

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
        fprintf(stdout, "Usage: AggregationTest [options]\n");
        fprintf(stdout, " --help             print this message\n");
        fprintf(stdout, " --verbose          run in verbose mode\n");
        fprintf(stdout, " --debug-level N    set lalDebugLevel\n");
        exit(0);
        break;

      case 'd':
        /* set debug level */
        lalDebugLevel = (INT4)optarg;
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

  return;
}


/* main program entry point */
INT4 main(INT4 argc, CHAR *argv[])
{
  /* declare variables */
  LIGOTimeGPS *start;
  LIGOTimeGPS *latest;
  CHAR *type;
  REAL8TimeSeries *series;
  INT4TimeSeries *dq_vector;
  INT4TimeSeries *state_vector;

  /* parameters */
  CHAR ifo[LIGOMETA_IFO_MAX] = "H1";
  LIGOTimeGPS gps = {918073010, 0};
  INT4 duration = 32;
  INT4 dq_bitmask = LAL_DQ_INJECTION;

  /* parse command line options */
  parse_options(argc, argv);

  /* determine gps time of latest frame file written */
  latest = XLALAggregationLatestGPS(ifo);
  if (latest == NULL)
  {
    static const char *func = "XLALAggregationLatestGPS";
    XLAL_ERROR(func, XLAL_EIO);
  }

  /* determine frame start time */
  start = XLALAggregationFrameStart(&gps);
  if (start == NULL)
  {
    static const char *func = "XLALAggregationFrameStart";
    XLAL_ERROR(func, XLAL_EINVAL);
  }
  /* check for correct value */
  if (start->gpsSeconds != 918073008)
  {
    static const char *func = "XLALAggregationFrameStart";
    XLAL_ERROR(func, XLAL_ETIME);
  }
  /* clear memory */
  XLALFree(start);

  /* determine frame type */
  type = XLALAggregationFrameType(ifo);
  if (type == NULL)
  {
    static const char *func = "XLALAggregationFrameType";
    XLAL_ERROR(func, XLAL_EINVAL);
  }
  /* check for correct value */
  if (strncmp(type, "H1_DMT_C00_L2", LIGOMETA_TYPE_MAX) != 0)
  {
    static const char *func = "XLALAggregationFrameType";
    XLAL_ERROR(func, XLAL_ENAME);
  }

  /* get strain data time series */
  series = XLALAggregationStrainData(ifo, &gps, duration);
  if (series == NULL)
  {
    fprintf(stderr, "failed: %d\n", xlalErrno);
    exit(xlalErrno);
  }
  LALDPrintTimeSeries(series, "series.dat");
  XLALDestroyREAL8TimeSeries(series);

  /* get data quality vector */
  dq_vector = XLALAggregationDQVector(ifo, &gps, duration);
  if (dq_vector == NULL)
  {
    fprintf(stderr, "failed: %d\n", xlalErrno);
    exit(xlalErrno);
  }
  LALI4PrintTimeSeries(dq_vector, "dq_vector.dat");
  XLALDestroyINT4TimeSeries(dq_vector);

  /* get state vector */
  state_vector = XLALAggregationStateVector(ifo, &gps, duration);
  if (state_vector == NULL)
  {
    fprintf(stderr, "failed: %d\n", xlalErrno);
    exit(xlalErrno);
  }
  LALI4PrintTimeSeries(state_vector, "state_vector.dat");
  XLALDestroyINT4TimeSeries(state_vector);

  /* get strain data time series, check data quality */
  series = XLALAggregationDQStrainData(ifo, &gps, duration, dq_bitmask);
  if (series == NULL)
  {
    fprintf(stderr, "failed: %d\n", xlalErrno);
    exit(xlalErrno);
  }
  LALDPrintTimeSeries(series, "series_dq.dat");
  XLALDestroyREAL8TimeSeries(series);

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  return 0;
}
