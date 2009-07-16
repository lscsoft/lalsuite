/*
 * fr_ninja.c - save numerical relativity waveforms as a frame
 *
 * Copyright (C) 2007, 2008 Adam Mercer
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
#include <lal/FrameStream.h>
#include <lal/LALDetectors.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include <FrameL.h>

#include <lalapps.h>

/* cvs info */
RCSID("$Id$");
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "fr_ninja"

/* defines */
/* TODO: how long can a FrHistory comment string be? */
#define HISTORY_COMMENT 512

/* minimum/maximum l values */
#define MIN_L 2
#define MAX_L 8

/* function prototypes */
static void print_usage(FILE *ptr, CHAR *program);


/* verbose flag */
extern int vrbflg;

/*
 * main program entry point
 */
INT4 main(INT4 argc, CHAR **argv)
{
  /* status */
  LALStatus status = blank_status;

  /* counters */
  int c;
  UINT4 i;
  
  /* mode counters */
  UINT4 l, m;

  /* metadata file/directory */
  CHAR *nrMetaFile = NULL;
  CHAR *nrDataDir = NULL;
  CHAR file_path[FILENAME_MAX];

  /* metadata parsing variables */
  LALParsedDataFile *meta_file = NULL;
  BOOLEAN wasRead = 0;
  CHAR *simulation_details = NULL;
  CHAR *nr_group = NULL;
  CHAR *email = NULL;
  CHAR *mass_ratio = NULL;
  CHAR *spin1x = NULL;
  CHAR *spin1y = NULL;
  CHAR *spin1z = NULL;
  CHAR *spin2x = NULL;
  CHAR *spin2y = NULL;
  CHAR *spin2z = NULL;
  CHAR *freqStart22 = NULL;
  CHAR *wf_name[MAX_L][(2*MAX_L) + 1];

  /* metadata */
  CHAR field[HISTORY_COMMENT];
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
  CHAR freq[HISTORY_COMMENT];
  CHAR creator[HISTORY_COMMENT];

  /* channel names */
  CHAR *plus_channel[MAX_L][(2*MAX_L) + 1];
  CHAR *cross_channel[MAX_L][(2*MAX_L) + 1];

  /* waveforms */
  UINT4 wf_length;
  REAL4TimeVectorSeries *waveforms[MAX_L][(2*MAX_L) + 1];
  REAL4TimeSeries *hplus[MAX_L][(2*MAX_L) + 1];
  REAL4TimeSeries *hcross[MAX_L][(2*MAX_L) + 1];

  /* frame variables */
  FrameH *frame;
  CHAR *frame_name = NULL;
  LIGOTimeGPS epoch;
  INT4 duration;
  INT4 detector_flags;

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

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level("33");

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
	fprintf( stdout, lalappsGitID );
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
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &freqStart22, meta_file, "freqStart22", &wasRead), &status);

  /* set waveform metadata */
  snprintf(sim, HISTORY_COMMENT, "simulation-details:%s", simulation_details);
  snprintf(group, HISTORY_COMMENT, "nr-group:%s", nr_group);
  snprintf(mail, HISTORY_COMMENT, "email:%s", email);
  snprintf(ratio, HISTORY_COMMENT, "mass-ratio:%s", mass_ratio);
  snprintf(s1x, HISTORY_COMMENT, "spin1x:%s", spin1x);
  snprintf(s1y, HISTORY_COMMENT, "spin1y:%s", spin1y);
  snprintf(s1z, HISTORY_COMMENT, "spin1z:%s", spin1z);
  snprintf(s2x, HISTORY_COMMENT, "spin2x:%s", spin2x);
  snprintf(s2y, HISTORY_COMMENT, "spin2y:%s", spin2y);
  snprintf(s2z, HISTORY_COMMENT, "spin2z:%s", spin2z);
  snprintf(freq, HISTORY_COMMENT, "freqStart22:%s", freqStart22);
  snprintf(creator, HISTORY_COMMENT, "creator:$Id$");

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
  XLALFrHistoryAdd(frame, "freqStart22", freq);
  XLALFrHistoryAdd(frame, "creator", creator);

  /* loop over l & m values */
  for (l = MIN_L; l <= MAX_L; l++)
  {
    for (m = (MAX_L - l); m <= MAX_L + l; m++)
    {
      /* ensure pointers are NULL */
      wf_name[l][m] = NULL;
      waveforms[l][m] = NULL;
      plus_channel[l][m] = NULL;
      cross_channel[l][m] = NULL;
      hplus[l][m] = NULL;
      hcross[l][m] = NULL;

      /* generate channel names */
      plus_channel[l][m] = XLALGetNinjaChannelName( "plus", l, m - MAX_L);
      cross_channel[l][m] = XLALGetNinjaChannelName("cross", l, m - MAX_L);

      /* initilise waveform time series */
      hplus[l][m] = XLALCreateREAL4TimeSeries(plus_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);
      hcross[l][m] = XLALCreateREAL4TimeSeries(cross_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);

      /* read ht-data section of metadata file */
      snprintf(field, HISTORY_COMMENT, "%d,%d", l, m - MAX_L);
      LAL_CALL(LALReadConfigSTRINGVariable(&status, &wf_name[l][m], meta_file, field, &wasRead), &status);

      /* read waveform */
      if (wf_name[l][m] != NULL)
      {
        /* get full path to waveform data file */
        snprintf(file_path, FILENAME_MAX, "%s/%s", nrDataDir, wf_name[l][m]);

        if (vrbflg)
          fprintf(stdout, "reading waveform: %s\n", file_path);

        /* read waveforms */
        LAL_CALL(LALReadNRWave_raw(&status, &waveforms[l][m], file_path), &status);
      }

      /* generate waveform time series from vector series */
      /* TODO: should use pointer arithmetic here and update the data
       * pointer in the REAL4TimeSeries to point to the appropriate
       * location within the REAL4TimeVector Series */
      if (waveforms[l][m])
      {
        /* get length of waveform */
        wf_length = waveforms[l][m]->data->vectorLength;

        /* allocate memory for waveform */
        XLALResizeREAL4TimeSeries(hplus[l][m], 0, wf_length);
        XLALResizeREAL4TimeSeries(hcross[l][m], 0, wf_length);

        /* set time spacing */
        hplus[l][m]->deltaT = waveforms[l][m]->deltaT;
        hcross[l][m]->deltaT = waveforms[l][m]->deltaT;

        /* copy waveforms into appropriate series */
        for (i = 0; i < wf_length; i ++) {
          hplus[l][m]->data->data[i] = waveforms[l][m]->data->data[i];
          hcross[l][m]->data->data[i] = waveforms[l][m]->data->data[wf_length + i];
        }
      }

     /* add channels to frame */
     if ((hplus[l][m]->data->length) && (hcross[l][m]->data->length))
      {
        XLALFrameAddREAL4TimeSeriesSimData(frame, hplus[l][m]);
        XLALFrameAddREAL4TimeSeriesSimData(frame, hcross[l][m]);
      }
    }
  }

  if (vrbflg)
    fprintf(stdout, "writing frame: %s\n", frame_name);

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
  for (l = MIN_L; l <= MAX_L; l++)
  {
    for (m = (MAX_L - l); m <= MAX_L + l; m++)
    {
      /* channel names */
      if (plus_channel[l][m])
        LALFree(plus_channel[l][m]);

      if (cross_channel[l][m])
        LALFree(cross_channel[l][m]);

      if (wf_name[l][m])
        LALFree(wf_name[l][m]);

      /* raw waveforms */
      if (waveforms[l][m]) {
        LALFree(waveforms[l][m]->data->data);
        LALFree(waveforms[l][m]->data);
        LALFree(waveforms[l][m]);
      }

      /* hplus */
      if (hplus[l][m])
        XLALDestroyREAL4TimeSeries(hplus[l][m]);

      /* hcross */
      if (hcross[l][m])
        XLALDestroyREAL4TimeSeries(hcross[l][m]);
    }
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
