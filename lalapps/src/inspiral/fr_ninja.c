/*
 * fr_ninja.c - save numerical relativity waveforms as a frame
 *
 * Copyright (C) 2007 Adam Mercer
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
 * Revision: $Id$
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>
#include <lal/LALFrameIO.h>
#include <lal/LALDetectors.h>
#include <lal/NRWaveIO.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lalapps.h>

/* cvs info */
RCSID("$Id$");
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "fr_ninja"

/* true/false */
#define TRUE {1==1}
#define FALSE {1==0}

/* defines */
/* TODO: how long can a FrHistory comment string be? */
#define HISTORY_COMMENT 512

/* function prototypes */
static void print_usage(FILE *ptr, CHAR *program);
static CHAR* channel_name(CHAR *polarisation, INT4 index, CHAR *channel);

/* verbose flag */
extern int vrbflg;

/* waveform enum */
enum {
  P2P2,
  P2P1,
  P2P0,
  P2N1,
  P2N2,
  NUM_WAVEFORMS};

/*
 * main program entry point
 */
INT4 main(INT4 argc, CHAR **argv)
{
  /* status */
  LALStatus status = blank_status;

  /* counters */
  int c;
  UINT4 i, j;

  /* metadata file/directory */
  CHAR *nrMetaFile = NULL;
  CHAR *nrDataDir = NULL;
  CHAR file_path[FILENAME_MAX];

  /* metadata parsing variables */
  LALParsedDataFile *meta_file = NULL;
  BOOLEAN wasRead = FALSE;
  CHAR *simulation_details = NULL;
  CHAR *nr_group = NULL;
  CHAR *email = NULL;
  CHAR *mass_ratio = NULL;;
  CHAR *spin1x = NULL;
  CHAR *spin1y = NULL;
  CHAR *spin1z = NULL;
  CHAR *spin2x = NULL;
  CHAR *spin2y = NULL;
  CHAR *spin2z = NULL;
  CHAR *wf_name[NUM_WAVEFORMS];

  /* metadata */
  CHAR sim[HISTORY_COMMENT];
  CHAR group[HISTORY_COMMENT];
  CHAR mail[HISTORY_COMMENT];
  CHAR ratio[HISTORY_COMMENT];
  CHAR s1x[HISTORY_COMMENT];
  CHAR s1y[HISTORY_COMMENT];
  CHAR s1z[HISTORY_COMMENT];
  CHAR s2x[HISTORY_COMMENT];
  CHAR s2y[HISTORY_COMMENT];
  CHAR s2z[HISTORY_COMMENT];

  /* channel names */
  CHAR *plus_channel[NUM_WAVEFORMS];
  CHAR *cross_channel[NUM_WAVEFORMS];

  /* waveforms */
  UINT4 wf_length;
  REAL4TimeVectorSeries *waveforms[NUM_WAVEFORMS];
  REAL4TimeSeries *hplus[NUM_WAVEFORMS];
  REAL4TimeSeries *hcross[NUM_WAVEFORMS];

  /* frame variables */
  FrameH *frame;
  CHAR *frame_name = NULL;
  LIGOTimeGPS epoch;
  INT4 duration;
  INT4 detector_flags;

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level("33");

  /* getopt arguments */
  struct option long_options[] =
  {
    /* options that set a flag */
    {"verbose", no_argument, &vrbflg, 1},
    /* options that don't set a flag */
    {"nr-meta-file", required_argument, 0, 'm'},
    {"nr-data-dir", required_argument, 0, 'd'},
    {"output", required_argument, 0, 'o'},
    {"debug-level", required_argument, 0, 'D'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'V'},
    {0, 0, 0, 0}
  };

  /* parse the arguments */
  while(1)
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, "m:d:o:D:hV", long_options, &option_index);

    /* detect the end of the options */
    if (c == -1)
      break;

    switch(c)
    {
      case 0:
        /* if this option sets a flag, do nothing else for now */
        if (long_options[option_index].flag != 0)
          break;
        else
        {
          fprintf(stderr, "Error parsing option %s with argument %s\n", long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'h':
        /* help message */
        print_usage(stdout, argv[0]);
        exit(0);
        break;

      case 'V':
        /* print version information and exit */
        fprintf(stdout, "Numerical Relativity Frame Generation\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n");
        exit(0);
        break;

      case 'm':
        /* create storage for the meta file name */
        optarg_len = strlen(optarg) + 1;
        nrMetaFile = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(nrMetaFile, optarg, optarg_len);
        break;

      case 'd':
        /* create storage for the meta data directory name */
        optarg_len = strlen(optarg) + 1;
        nrDataDir = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(nrDataDir, optarg, optarg_len);
        break;

      case 'o':
        /* create storage for the output frame file name */
        optarg_len = strlen(optarg) + 1;
        frame_name = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(frame_name, optarg, optarg_len);
        break;

      case 'D':
        /* set debug level */
        set_debug_level(optarg);
        break;

      case '?':
        print_usage(stderr, argv[0]);
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing arguments\n");
        print_usage(stderr, argv[0]);
        exit(1);
        break;
    }
  }

  /* check for extraneous command line arguments */
  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
      fprintf(stderr, "%s\n", argv[optind++]);
    exit(1);
  }

  /*
   * check validity of arguments
   */

  /* meta file specified */
  if (nrMetaFile == NULL)
  {
    fprintf(stderr, "--nr-meta-file must be specified\n");
    exit(1);
  }

  /* data directory specified */
  if (nrDataDir == NULL)
  {
    fprintf(stderr, "--nr-data-dir must be specified\n");
    exit(1);
  }

  /* output frame filename specified */
  if (frame_name == NULL)
  {
    fprintf(stderr, "--output must be specified\n");
    exit(1);
  }

  /*
   * main code
   */

  /* frame metadata */
  /* TODO: set these to something sensible */
  duration = 0;
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;
  detector_flags = 0;

  /* initialise waveform pointers */
  for (i = 0; i < NUM_WAVEFORMS; i++)
  {
    /* ensure pointers are NULL */
    wf_name[i] = NULL;
    waveforms[i] = NULL;
    plus_channel[i] = NULL;
    cross_channel[i] = NULL;
    hplus[i] = NULL;
    hcross[i] = NULL;

    /* generate channel names */
    plus_channel[i] = channel_name("plus", i, plus_channel[i]);
    cross_channel[i] = channel_name("cross", i, cross_channel[i]);

    /* initilise waveform time series */
    hplus[i] = XLALCreateREAL4TimeSeries(plus_channel[i], &epoch, 0, 0, &lalDimensionlessUnit, 0);
    hcross[i] = XLALCreateREAL4TimeSeries(cross_channel[i], &epoch, 0, 0, &lalDimensionlessUnit, 0);
  }

  if (vrbflg)
    fprintf(stdout, "reading metadata: %s\n", nrMetaFile);

  /* open metadata file */
  LAL_CALL(LALParseDataFile(&status, &meta_file, nrMetaFile), &status);

  /* metadata section */
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &simulation_details, meta_file, "simulation-details", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &nr_group, meta_file, "nr-group", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &email, meta_file, "email", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &mass_ratio, meta_file, "mass-ratio", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin1x, meta_file, "spin1x", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin1y, meta_file, "spin1y", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin1z, meta_file, "spin1z", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin2x, meta_file, "spin2x", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin2y, meta_file, "spin2y", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &spin2z, meta_file, "spin2z", &wasRead), &status);

  /* ht-data section */
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[P2P2], meta_file, "2,2", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[P2P1], meta_file, "2,1", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[P2P0], meta_file, "2,0", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[P2N1], meta_file, "2,-1", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[P2N2], meta_file, "2,-2", &wasRead), &status);

  /* read waveforms */
  for (i = 0; i < NUM_WAVEFORMS; i++)
  {
    /* get full path to waveform data file */
    LALSnprintf(file_path, FILENAME_MAX, "%s/%s", nrDataDir, wf_name[i]);

    if (vrbflg)
      fprintf(stdout, "reading waveform: %s\n", file_path);

    /* read waveforms */
    LAL_CALL(LALReadNRWave_raw(&status, &waveforms[i], file_path), &status);
  }

  /* set waveform metadata */
  LALSnprintf(sim, HISTORY_COMMENT, "simulation-details:%s", simulation_details);
  LALSnprintf(group, HISTORY_COMMENT, "nr-group:%s", nr_group);
  LALSnprintf(mail, HISTORY_COMMENT, "email:%s", email);
  LALSnprintf(ratio, HISTORY_COMMENT, "mass-ratio:%s", mass_ratio);
  LALSnprintf(s1x, HISTORY_COMMENT, "spin1x:%s", spin1x);
  LALSnprintf(s1y, HISTORY_COMMENT, "spin1y:%s", spin1y);
  LALSnprintf(s1z, HISTORY_COMMENT, "spin1z:%s", spin1z);
  LALSnprintf(s2x, HISTORY_COMMENT, "spin2x:%s", spin2x);
  LALSnprintf(s2y, HISTORY_COMMENT, "spin2y:%s", spin2y);
  LALSnprintf(s2z, HISTORY_COMMENT, "spin2z:%s", spin2z);

  /* generate waveform time series from vector series */
  /* TODO: should use pointer arithmetic here and update the data
   * pointer in the REAL4TimeSeries to point to the appropriate
   * location within the REAL4TimeVector Series */
  for (i = 0; i < NUM_WAVEFORMS; i++)
  {
    /* get length of waveform */
    wf_length = waveforms[i]->data->vectorLength;

    /* allocate memory for waveform */
    XLALResizeREAL4TimeSeries(hplus[i], 0, wf_length);
    XLALResizeREAL4TimeSeries(hcross[i], 0, wf_length);

    /* set time spacing */
    hplus[i]->deltaT = waveforms[i]->deltaT;
    hcross[i]->deltaT = waveforms[i]->deltaT;

    /* copy waveforms into appropriate series */
    for (j = 0; j < wf_length; j ++)
    {
      hplus[i]->data->data[j] = waveforms[i]->data->data[j];
      hcross[i]->data->data[j] = waveforms[i]->data->data[wf_length + j];
    }
  }

  if (vrbflg)
    fprintf(stdout, "creating and saving frame: %s\n", frame_name);

  /* define frame */
  frame = XLALFrameNew(&epoch, duration, "NR", 0, 1, detector_flags);

  /* add metadata as FrHistory structures */
  XLALFrHistoryAdd(frame, "simulation-details", sim);
  XLALFrHistoryAdd(frame, "nr-group", group);
  XLALFrHistoryAdd(frame, "email", mail);
  XLALFrHistoryAdd(frame, "mass-ratio", ratio);
  XLALFrHistoryAdd(frame, "spin1x", s1x);
  XLALFrHistoryAdd(frame, "spin1y", s1y);
  XLALFrHistoryAdd(frame, "spin1z", s1z);
  XLALFrHistoryAdd(frame, "spin2x", s2x);
  XLALFrHistoryAdd(frame, "spin2y", s2y);
  XLALFrHistoryAdd(frame, "spin2z", s2z);

  /* add channels to frame */
  for (i = 0; i < NUM_WAVEFORMS; i++)
  {
    XLALFrameAddREAL4TimeSeriesSimData(frame, hplus[i]);
    XLALFrameAddREAL4TimeSeriesSimData(frame, hcross[i]);
  }

  /* write frame */
  if (XLALFrameWrite(frame, frame_name, 8) != 0 )
  {
    fprintf(stderr, "Error: Cannot save frame file '%s'\n", frame_name);
    exit(1);
  }

  /*
   * clear memory
   */

  /* strings */
  free(nrMetaFile);
  free(nrDataDir);
  free(frame_name);
  LALFree(simulation_details);
  LALFree(nr_group);
  LALFree(email);
  LALFree(mass_ratio);
  LALFree(spin1x);
  LALFree(spin1y);
  LALFree(spin1z);
  LALFree(spin2x);
  LALFree(spin2y);
  LALFree(spin2z);

  /* config file */
  LALFree(meta_file->lines->list->data);
  LALFree(meta_file->lines->list);
  LALFree(meta_file->lines->tokens);
  LALFree(meta_file->lines);
  LALFree(meta_file->wasRead);
  LALFree(meta_file);

  /* waveforms */
  for (i = 0; i < NUM_WAVEFORMS; i++)
  {
    /* channel names */
    LALFree(plus_channel[i]);
    LALFree(cross_channel[i]);
    LALFree(wf_name[i]);

    /* raw waveforms */
    LALFree(waveforms[i]->data->data);
    LALFree(waveforms[i]->data);
    LALFree(waveforms[i]);

    /* hplus */
    LALFree(hplus[i]->data->data);
    LALFree(hplus[i]->data);
    LALFree(hplus[i]);

    /* hcross */
    LALFree(hcross[i]->data->data);
    LALFree(hcross[i]->data);
    LALFree(hcross[i]);
  }

  /* clear frame */
  FrameFree(frame);

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  exit(0);
}


/* function to display usage information */
static void print_usage(FILE *ptr, CHAR *program)
{
  fprintf(ptr,
      "Usage: %s [options\n"\
      "[--help                    display this message and exit]\n"\
      "[--version                 print version information and exit]\n"\
      "[--verbose                 display progress information]\n"\
      "[--debug-level    LEVEL    set the debug level]\n"\
      " --nr-meta-file   FILE     file containing the details of the available\n"\
      "                           numerical relativity waveforms\n"\
      " --nr-data-dir    DIR      directory containing the numerical relativity\n"\
      "                           waveforms\n"\
      " --output         FILE     name of output frame file\n", program);
}


/* function to return channel name */
static CHAR* channel_name(CHAR *polarisation, INT4 index, CHAR *channel)
{
  /* check that channel is a NULL pointer */
  if (channel != NULL)
  {
    fprintf(stderr, "Error: 'channel' should be a NULL pointer\n");
    exit(1);
  }

  /* check for suppored polarisation */
  if (strncmp(polarisation, "plus", 4) == 0)
  {
    /* supported polarisation */
  }
  else if (strncmp(polarisation, "cross", 5) == 0)
  {
    /* supported polarisation */
  }
  else
  {
    fprintf(stderr, "Error: Unknown polarisation '%s'\n", polarisation);
    exit(1);
  }

  /* allocate memory for channel */
  channel = (CHAR *)LALCalloc(LIGOMETA_CHANNEL_MAX, sizeof(CHAR));

  /* construct channel name */
  if (index == P2P2)
    LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l2_mp2", polarisation);
  else if (index == P2P1)
    LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l2_mp1", polarisation);
  else if (index == P2P0)
    LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l2_mp0", polarisation);
  else if (index == P2N1)
    LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l2_mn1", polarisation);
  else if (index == P2N2)
    LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l2_mn2", polarisation);
  else
  {
    fprintf(stderr, "Error: Unknown waveform index: '%d'\n", index);
    exit(1);
  }

  /* return channel name */
  return channel;
} 
