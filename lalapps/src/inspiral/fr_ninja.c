/*
 * fr_ninja.c - save numerical relativity waveforms as a frame
 *
 * Copyright (C) 2007,2008,2010 Adam Mercer
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
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>
#include <lal/LALDetectors.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALFrameL.h>

#include <lalapps.h>
#include <LALAppsVCSInfo.h>

/* program info */
#define PROGRAM_NAME "lalapps_fr_ninja"

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

  /* metadata format */
  CHAR *metadata_format = NULL;

  /* metadata parsing variables */
  LALParsedDataFile *meta_file = NULL;
  BOOLEAN wasRead = 0;
  CHAR field[HISTORY_COMMENT];
  CHAR *wf_name[MAX_L+1][(2*MAX_L) + 1];

  /* common metadata */
  CHAR *md_mass_ratio = NULL;
  CHAR *md_spin1x = NULL;
  CHAR *md_spin1y = NULL;
  CHAR *md_spin1z = NULL;
  CHAR *md_spin2x = NULL;
  CHAR *md_spin2y = NULL;
  CHAR *md_spin2z = NULL;
  CHAR *md_freq_start_22 = NULL;

  /* NINJA1 metadata */
  CHAR *md_simulation_details = NULL;
  CHAR *md_nr_group = NULL;
  CHAR *md_email = NULL;

  /* NINJA2 metadata */
  CHAR *md_waveform_name = NULL;
  CHAR *md_initial_separation = NULL;
  CHAR *md_eccentricity = NULL;
  CHAR *md_number_of_cycles_22 = NULL;
  CHAR *md_code = NULL;
  CHAR *md_submitter_email = NULL;
  CHAR *md_authors_emails = NULL;

  /* common metadata strings */
  CHAR str_mass_ratio[HISTORY_COMMENT];
  CHAR str_spin1x[HISTORY_COMMENT];
  CHAR str_spin1y[HISTORY_COMMENT];
  CHAR str_spin1z[HISTORY_COMMENT];
  CHAR str_spin2x[HISTORY_COMMENT];
  CHAR str_spin2y[HISTORY_COMMENT];
  CHAR str_spin2z[HISTORY_COMMENT];
  CHAR str_freq_start_22[HISTORY_COMMENT];
  CHAR str_creator[HISTORY_COMMENT];

  /* NINJA1 metadata strings */
  CHAR str_simulation_details[HISTORY_COMMENT];
  CHAR str_nr_group[HISTORY_COMMENT];
  CHAR str_email[HISTORY_COMMENT];

  /* NINJA2 metadata strings */
  CHAR str_waveform_name[HISTORY_COMMENT];
  CHAR str_initial_separation[HISTORY_COMMENT];
  CHAR str_eccentricity[HISTORY_COMMENT];
  CHAR str_number_of_cycles_22[HISTORY_COMMENT];
  CHAR str_code[HISTORY_COMMENT];
  CHAR str_submitter_email[HISTORY_COMMENT];
  CHAR str_authors_emails[HISTORY_COMMENT];

  /* channel names */
  CHAR *plus_channel[MAX_L+1][(2*MAX_L) + 1];
  CHAR *cross_channel[MAX_L+1][(2*MAX_L) + 1];

  /* waveforms */
  UINT4 wf_length;
  REAL4TimeVectorSeries *waveforms[MAX_L][(2*MAX_L) + 1];
  REAL4TimeSeries *hplus[MAX_L+1][(2*MAX_L) + 1];
  REAL4TimeSeries *hcross[MAX_L+1][(2*MAX_L) + 1];

  REAL8TimeVectorSeries *waveformsREAL8[MAX_L][(2*MAX_L) + 1];
  REAL8TimeSeries *hplusREAL8[MAX_L+1][(2*MAX_L) + 1];
  REAL8TimeSeries *hcrossREAL8[MAX_L+1][(2*MAX_L) + 1];

  /* frame variables */
  LALFrameH *frame;
  CHAR *frame_name = NULL;
  LIGOTimeGPS epoch;
  INT4 duration;
  INT4 detector_flags;

  INT4 generatingREAL8 = 0;

  /* getopt arguments */
  struct option long_options[] =
  {
    /* options that set a flag */
    {"verbose", no_argument, &vrbflg, 1},
    /* options that don't set a flag */
    {"format", required_argument, 0, 'f'},
    {"nr-meta-file", required_argument, 0, 'm'},
    {"nr-data-dir", required_argument, 0, 'd'},
    {"output", required_argument, 0, 'o'},
    {"double-precision", no_argument, 0, 'p'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'V'},
    {0, 0, 0, 0}
  };

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;

  /* parse the arguments */
  while(1)
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, "f:m:d:o:phV", long_options, &option_index);

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
        fprintf(stdout, "Numerical Relativity Frame Generation\n");
        XLALOutputVersionString(stderr, 0);
        exit(0);
        break;

      case 'f':
        /* create storage for the metadata format */
        optarg_len = strlen(optarg) + 1;
        metadata_format = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(metadata_format, optarg, optarg_len);
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

      case 'p':
        /* We're generating a double-precision frame */
        generatingREAL8 = 1;
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

  /* metadata format specified */
  if (metadata_format == NULL)
  {
    fprintf(stderr, "warning: --format not specified, assuming NINJA1\n");
    metadata_format = (CHAR *)calloc(7, sizeof(CHAR));
    memcpy(metadata_format, "NINJA1", 7);
  }

  /* check for supported metadata format */
  if (strcmp(metadata_format, "NINJA1") == 0);
  else if (strcmp(metadata_format, "NINJA2") == 0);
  else
  {
    fprintf(stderr, "Supported metadata formats are NINJA1 and NINJA2 (%s specified)\n", metadata_format);
    exit(1);
  }

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

  /*
   * get metadata
   */

  /* common metadata */
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_mass_ratio, meta_file, "mass-ratio", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin1x, meta_file, "spin1x", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin1y, meta_file, "spin1y", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin1z, meta_file, "spin1z", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin2x, meta_file, "spin2x", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin2y, meta_file, "spin2y", &wasRead), &status);
  LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_spin2z, meta_file, "spin2z", &wasRead), &status);

  /* format specific metadata */
  if (strcmp(metadata_format, "NINJA1") == 0)
  {
    /* NINJA1 */
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_simulation_details, meta_file, "simulation-details", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_nr_group, meta_file, "nr-group", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_email, meta_file, "email", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_freq_start_22, meta_file, "freqStart22", &wasRead), &status);
  }
  else if (strcmp(metadata_format, "NINJA2") == 0)
  {
    /* NINJA2 */
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_waveform_name, meta_file, "waveform-name", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_initial_separation, meta_file, "initial-separation", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_eccentricity, meta_file, "eccentricity", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_number_of_cycles_22, meta_file, "number-of-cycles-22", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_code, meta_file, "code", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_submitter_email, meta_file, "submitter-email", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_authors_emails, meta_file, "authors-emails", &wasRead), &status);
    LAL_CALL(LALReadConfigSTRINGVariable(&status, &md_freq_start_22, meta_file, "freq-start-22", &wasRead), &status);
  }
  else
  {
    /* unknown metadata format - should not be executed */
    fprintf(stderr, "error: unsupported metadata format: %s\n", metadata_format);
    exit(1);
  }

  /*
   * set metadata strings
   */

  /* common waveform */
  snprintf(str_mass_ratio, HISTORY_COMMENT, "mass-ratio:%s", md_mass_ratio);
  snprintf(str_spin1x, HISTORY_COMMENT, "spin1x:%s", md_spin1x);
  snprintf(str_spin1y, HISTORY_COMMENT, "spin1y:%s", md_spin1y);
  snprintf(str_spin1z, HISTORY_COMMENT, "spin1z:%s", md_spin1z);
  snprintf(str_spin2x, HISTORY_COMMENT, "spin2x:%s", md_spin2x);
  snprintf(str_spin2y, HISTORY_COMMENT, "spin2y:%s", md_spin2y);
  snprintf(str_spin2z, HISTORY_COMMENT, "spin2z:%s", md_spin2z);
  snprintf(str_creator, HISTORY_COMMENT, "creator:%s(git:%s)", PROGRAM_NAME, lalAppsVCSId);

  /* format specific metadata */
  if (strcmp(metadata_format, "NINJA1") == 0)
  {
    /* NINJA1 */
    snprintf(str_freq_start_22, HISTORY_COMMENT, "freqStart22:%s", md_freq_start_22);
    snprintf(str_simulation_details, HISTORY_COMMENT, "simulation-details:%s", md_simulation_details);
    snprintf(str_nr_group, HISTORY_COMMENT, "nr-group:%s", md_nr_group);
    snprintf(str_email, HISTORY_COMMENT, "email:%s", md_email);
  }
  else if (strcmp(metadata_format, "NINJA2") == 0)
  {
    /* NINJA2 */
    snprintf(str_waveform_name, HISTORY_COMMENT, "waveform-name:%s", md_waveform_name);
    snprintf(str_initial_separation, HISTORY_COMMENT, "inital-separation:%s", md_initial_separation);
    snprintf(str_eccentricity, HISTORY_COMMENT, "eccentricity:%s", md_eccentricity);
    snprintf(str_freq_start_22, HISTORY_COMMENT, "freq_start_22:%s", md_freq_start_22);
    snprintf(str_number_of_cycles_22, HISTORY_COMMENT, "number-of-cycles-22:%s", md_number_of_cycles_22);
    snprintf(str_code, HISTORY_COMMENT, "code:%s", md_code);
    snprintf(str_submitter_email, HISTORY_COMMENT, "submitter-email:%s", md_submitter_email);
    snprintf(str_authors_emails, HISTORY_COMMENT, "authors-emails:%s", md_authors_emails);
  }
  else
  {
    /* unknown metadata format - should not be executed */
    fprintf(stderr, "error: unsupported metadata format: %s\n", metadata_format);
    exit(1);
  }

  /* define frame */
  frame = XLALFrameNew(&epoch, duration, "NR", 0, 1, detector_flags);

  /*
   * add metadata as FrHistory structures
   */

  /* common metadata */
  XLALFrameAddFrHistory(frame, "creator", str_creator);
  XLALFrameAddFrHistory(frame, "mass-ratio", str_mass_ratio);
  XLALFrameAddFrHistory(frame, "spin1x", str_spin1x);
  XLALFrameAddFrHistory(frame, "spin1y", str_spin1y);
  XLALFrameAddFrHistory(frame, "spin1z", str_spin1z);
  XLALFrameAddFrHistory(frame, "spin2x", str_spin2x);
  XLALFrameAddFrHistory(frame, "spin2y", str_spin2y);
  XLALFrameAddFrHistory(frame, "spin2z", str_spin2z);

  /* format specific metadata */
  if (strcmp(metadata_format, "NINJA1") == 0)
  {
    /* NINJA1 */
    XLALFrameAddFrHistory(frame, "simulation-details", str_simulation_details);
    XLALFrameAddFrHistory(frame, "nr-group", str_nr_group);
    XLALFrameAddFrHistory(frame, "email", str_email);
    XLALFrameAddFrHistory(frame, "freqStart22", str_freq_start_22);
  }
  else if (strcmp(metadata_format, "NINJA2") == 0)
  {
    /* NINJA2 */
    XLALFrameAddFrHistory(frame, "waveform-name", str_waveform_name);
    XLALFrameAddFrHistory(frame, "initial-separation", str_initial_separation);
    XLALFrameAddFrHistory(frame, "eccentricity", str_eccentricity);
    XLALFrameAddFrHistory(frame, "freq_start_22", str_freq_start_22);
    XLALFrameAddFrHistory(frame, "number-of-cycles-22", str_number_of_cycles_22);
    XLALFrameAddFrHistory(frame, "code", str_code);
    XLALFrameAddFrHistory(frame, "submitter-email", str_code);
    XLALFrameAddFrHistory(frame, "authors-emails", str_authors_emails);
  }
  else
  {
    /* unknown metadata format - should not be executed */
    fprintf(stderr, "error: unsupported metadata format: %s\n", metadata_format);
    exit(1);
  }

  /* loop over l & m values */
  for (l = MIN_L; l <= MAX_L; l++)
  {
    for (m = (MAX_L - l); m <= MAX_L + l; m++)
    {
      /* ensure pointers are NULL */
      wf_name[l][m] = NULL;
      plus_channel[l][m] = NULL;
      cross_channel[l][m] = NULL;

      /* generate channel names */
      plus_channel[l][m] = XLALGetNinjaChannelName("plus", l, m - MAX_L);
      cross_channel[l][m] = XLALGetNinjaChannelName("cross", l, m - MAX_L);

      if (generatingREAL8)
      {
        hplusREAL8[l][m] = NULL;
        hcrossREAL8[l][m] = NULL;
        waveformsREAL8[l][m] = NULL;
        
        /* initilise waveform time series */
        hplusREAL8[l][m] = XLALCreateREAL8TimeSeries(plus_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);
        hcrossREAL8[l][m] = XLALCreateREAL8TimeSeries(cross_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);
      }
      else
      {
        hplus[l][m] = NULL;
        hcross[l][m] = NULL;
        waveforms[l][m] = NULL;
      
        /* initilise waveform time series */
        hplus[l][m] = XLALCreateREAL4TimeSeries(plus_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);
        hcross[l][m] = XLALCreateREAL4TimeSeries(cross_channel[l][m], &epoch, 0, 0, &lalDimensionlessUnit, 0);
      }

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
        if (generatingREAL8)
        {
          LAL_CALL(LALReadNRWave_raw_real8(&status, &waveformsREAL8[l][m], file_path), &status);
        }
        else
        {
          LAL_CALL(LALReadNRWave_raw(&status, &waveforms[l][m], file_path), &status);
        }
      }

      /* generate waveform time series from vector series */
      /* TODO: should use pointer arithmetic here and update the data
       * pointer in the REAL4TimeSeries to point to the appropriate
       * location within the REAL4TimeVector Series */
      if (generatingREAL8) {
        if (waveformsREAL8[l][m])
        {
          /* get length of waveform */
          wf_length = waveformsREAL8[l][m]->data->vectorLength;

          /* allocate memory for waveform */
          XLALResizeREAL8TimeSeries(hplusREAL8[l][m], 0, wf_length);
          XLALResizeREAL8TimeSeries(hcrossREAL8[l][m], 0, wf_length);

          /* set time spacing */
          hplusREAL8[l][m]->deltaT = waveformsREAL8[l][m]->deltaT;
          hcrossREAL8[l][m]->deltaT = waveformsREAL8[l][m]->deltaT;

          /* copy waveforms into appropriate series */
          for (i = 0; i < wf_length; i ++) {
            hplusREAL8[l][m]->data->data[i] = waveformsREAL8[l][m]->data->data[i];
            hcrossREAL8[l][m]->data->data[i] = waveformsREAL8[l][m]->data->data[wf_length + i];
          }

          /* Done with waveformsREAL8, clean up here to limit memory usage */
          LALFree(waveformsREAL8[l][m]->data->data);
          LALFree(waveformsREAL8[l][m]->data);
          LALFree(waveformsREAL8[l][m]);
        }
      }
      else /* REAL4 */
      {
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
      }

      /* add channels to frame */
      if (generatingREAL8)
      {
        if ((hplusREAL8[l][m]->data->length) && (hcrossREAL8[l][m]->data->length))
        {
          XLALFrameAddREAL8TimeSeriesSimData(frame, hplusREAL8[l][m]);
          XLALFrameAddREAL8TimeSeriesSimData(frame, hcrossREAL8[l][m]);
        }
      }
      else
      {
        if ((hplus[l][m]->data->length) && (hcross[l][m]->data->length))
        {
          XLALFrameAddREAL4TimeSeriesSimData(frame, hplus[l][m]);
          XLALFrameAddREAL4TimeSeriesSimData(frame, hcross[l][m]);
        }
      }
    }
  }

  if (vrbflg)
    fprintf(stdout, "writing frame: %s\n", frame_name);

  /* write frame */
  if (XLALFrameWrite(frame, frame_name) != 0 )
  {
    fprintf(stderr, "Error: Cannot save frame file '%s'\n", frame_name);
    exit(1);
  }

  /*
   * clear memory
   */

  /* strings */
  if(nrMetaFile) free(nrMetaFile);
  if(nrDataDir)  free(nrDataDir);
  if(frame_name) free(frame_name);
  if(metadata_format) free(metadata_format);

  /* common metadata */
  if(md_mass_ratio) LALFree(md_mass_ratio);
  if(md_spin1x) LALFree(md_spin1x);
  if(md_spin1y) LALFree(md_spin1y);
  if(md_spin1z) LALFree(md_spin1z);
  if(md_spin2x) LALFree(md_spin2x);
  if(md_spin2y) LALFree(md_spin2y);
  if(md_spin2z) LALFree(md_spin2z);
  if(md_freq_start_22) LALFree(md_freq_start_22);

  /* NINJA1 metadata */
  if(md_simulation_details) LALFree(md_simulation_details);
  if(md_nr_group)           LALFree(md_nr_group);
  if(md_email)              LALFree(md_email);

  /* NINJA2 metadata */
  if(md_waveform_name)       LALFree(md_waveform_name);
  if(md_initial_separation)  LALFree(md_initial_separation);
  if(md_eccentricity)        LALFree(md_eccentricity);
  if(md_number_of_cycles_22) LALFree(md_number_of_cycles_22);
  if(md_code)                LALFree(md_code);
  if(md_submitter_email)     LALFree(md_submitter_email);
  if(md_authors_emails)      LALFree(md_authors_emails);

  /* config file */
  if(meta_file->lines->list->data) LALFree(meta_file->lines->list->data);
  if(meta_file->lines->list)   LALFree(meta_file->lines->list);
  if(meta_file->lines->tokens) LALFree(meta_file->lines->tokens);
  if(meta_file->lines)   LALFree(meta_file->lines);
  if(meta_file->wasRead) LALFree(meta_file->wasRead);
  if(meta_file)          LALFree(meta_file);

  /* waveforms */
  if (generatingREAL8)
  {
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

        /* hplus */
        if (hplusREAL8[l][m])
          XLALDestroyREAL8TimeSeries(hplusREAL8[l][m]);

        /* hcross */
        if (hcrossREAL8[l][m])
          XLALDestroyREAL8TimeSeries(hcrossREAL8[l][m]);
      }
    }
  }
  else
  {
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
  }

  /* clear frame */
  XLALFrameFree(frame);

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
      "[--format         FORMAT   metadata format, defaults to NINJA1]\n"\
      "[--double-precision        generate REAL8 files, default is REAL4]\n"\
      " --nr-meta-file   FILE     file containing the details of the available\n"\
      "                           numerical relativity waveforms\n"\
      " --nr-data-dir    DIR      directory containing the numerical relativity\n"\
      "                           waveforms\n"\
      " --output         FILE     name of output frame file\n", program);
}
