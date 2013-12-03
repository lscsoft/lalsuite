/*
 * stochastic.c - SGWB Standalone Analysis Pipeline
 *
 * Copyright (C) 2002-2006,2009,2010 Adam Mercer
 * Copyright (C) 2003-2004 Tania Regimbau
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#include "data_input.h"
#include "misc.h"
#include "sgwb.h"
#include "data_output.h"
#include "stochastic.h"


/* program info */
const CHAR* prog_name="lalapps_stochastic";

/* flags for getopt_long */
int middle_segment_flag;
int apply_mask_flag;
int high_pass_flag;
int overlap_hann_flag;
int recentre_flag;
int cc_spectra_flag;
int debug_flag;
extern int vrbflg;

/* xml comment/tags */
CHAR comment[LIGOMETA_COMMENT_MAX];
CHAR *userTag = NULL;

/* xml tables */
MetadataTable proctable;
MetadataTable procparams;
ProcessParamsTable *this_proc_param;

/* parameters for the stochastic search */

/* sampling parameters */
INT4 resampleRate;

/* data parameters */
INT4 startTime = 0;
INT4 endTime = 0;
INT4 intervalDuration = -1;
INT4 segmentDuration = -1;
INT4 calibOffset = -1;
CHAR *frameCacheOne = NULL;
CHAR *frameCacheTwo = NULL;
CHAR *calCacheOne = NULL;
CHAR *calCacheTwo = NULL;
CHAR *channelOne = NULL;
CHAR *channelTwo = NULL;
CHAR *ifoOne = NULL;
CHAR *ifoTwo = NULL;
INT4 siteOne;
INT4 siteTwo;

/* frequency band */
REAL8 fMin = -1;
REAL8 fMax = -1;

/* omegaGW parameters */
REAL4 alpha = 0;
REAL4 fRef = 100;
REAL4 omegaRef = 1;

/* window parameters */
INT4 hannDuration = -1;

/* high pass filtering parameters */
REAL4 highPassFreq = -1;
REAL4 highPassAtten = -1;
INT4  highPassOrder = -1;

/* GEO scale factor */
REAL4 geoScaleFactor = 1e18;

/* GEO high pass filter parameters */
REAL4 geoHighPassFreq = -1;
INT4  geoHighPassOrder = -1;
REAL4 geoHighPassAtten = -1;

/* number of bins for frequency masking */
INT4 maskBin = -1;

/* output file */
CHAR *outputPath = NULL;

/* helper functions */

/* display usage information */
static void display_usage(void)
{
  fprintf(stdout, "Usage: %s [options]\n", prog_name);
  fprintf(stdout, " --help                        print this message\n");
  fprintf(stdout, " --version                     display version\n");
  fprintf(stdout, " --verbose                     verbose mode\n");
  fprintf(stdout, " --debug                       debug mode\n");
  fprintf(stdout, " --user-tag STRING             set the user tag\n");
  fprintf(stdout, " --comment STRING              set the comment\n");
  fprintf(stdout, " --output-dir DIR              directory for output files\n");
  fprintf(stdout, " --cc-spectra                  save out cross correlation spectra\n");
  fprintf(stdout, " --gps-start-time N            GPS start time\n");
  fprintf(stdout, " --gps-end-time N              GPS end time\n");
  fprintf(stdout, " --interval-duration N         interval duration\n");
  fprintf(stdout, " --segment-duration N          segment duration\n");
  fprintf(stdout, " --resample-rate N             resample rate\n");
  fprintf(stdout, " --f-min N                     minimal frequency\n");
  fprintf(stdout, " --f-max N                     maximal frequency\n");
  fprintf(stdout, " --ifo-one IFO                 ifo for first stream\n");
  fprintf(stdout, " --ifo-two IFO                 ifo for second stream\n");
  fprintf(stdout, " --channel-one CHANNEL         channel for first stream\n");
  fprintf(stdout, " --channel-two CHANNEL         channel for second stream\n");
  fprintf(stdout, " --frame-cache-one FILE        cache file for first stream\n");
  fprintf(stdout, " --frame-cache-two FILE        cache file for second stream\n");
  fprintf(stdout, " --calibration-cache-one FILE  first stream calibration cache\n");
  fprintf(stdout, " --calibration-cache-two FILE  second stream calibration cache\n");
  fprintf(stdout, " --calibration-offset N        calibration offset\n");
  fprintf(stdout, " --apply-mask                  apply frequency masking\n");
  fprintf(stdout, " --mask-bin N                  number of bins to mask\n");
  fprintf(stdout, " --overlap-hann                overlaping hann windows\n");
  fprintf(stdout, " --hann-duration N             hann duration\n");
  fprintf(stdout, " --high-pass-filter            apply high pass filtering\n");
  fprintf(stdout, " --hpf-frequency N             high pass filter knee frequency\n");
  fprintf(stdout, " --hpf-attenuation N           high pass filter attenuation\n");
  fprintf(stdout, " --hpf-order N                 high pass filter order\n");
  fprintf(stdout, " --recentre                    recentre jobs\n");
  fprintf(stdout, " --middle-segment              use middle segment in PSD estimation\n");
  fprintf(stdout, " --geo-hpf-frequency N         GEO high pass filter knee frequency\n");
  fprintf(stdout, " --geo-hpf-attenuation N       GEO high pass filter attenuation\n");
  fprintf(stdout, " --geo-hpf-order N             GEO high pass filter order\n");
  fprintf(stdout, " --alpha N                     exponent on filter spectrum\n");
  fprintf(stdout, " --f-ref N                     reference frequency for filter spectrum\n");
  fprintf(stdout, " --omega0 N                    reference omega_0 for filter spectrum\n");
}

/* parse command line options */
static void parse_options(INT4 argc, CHAR *argv[])
{
  int c = -1;
  struct stat fileStatus;

  /* tempory variables */
  CHAR *channelOneTemp = NULL;
  CHAR *channelTwoTemp = NULL;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"middle-segment", no_argument, &middle_segment_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &vrbflg, 1},
      {"recentre", no_argument, &recentre_flag, 1},
      {"cc-spectra", no_argument, &cc_spectra_flag, 1},
      {"debug", no_argument, &debug_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'a'},
      {"version", no_argument, 0, 'b'},
      {"user-tag", required_argument, 0, 'd'},
      {"comment", required_argument, 0, 'e'},
      {"output-dir", required_argument, 0, 'f'},
      {"gps-start-time", required_argument, 0, 'g'},
      {"gps-end-time", required_argument, 0, 'h'},
      {"interval-duration", required_argument, 0, 'i'},
      {"segment-duration", required_argument, 0, 'j'},
      {"resample-rate", required_argument, 0, 'k'},
      {"f-min", required_argument, 0, 'l'},
      {"f-max", required_argument, 0, 'm'},
      {"ifo-one", required_argument, 0, 'n'},
      {"ifo-two", required_argument, 0, 'o'},
      {"channel-one", required_argument, 0, 'p'},
      {"channel-two", required_argument, 0, 'q'},
      {"frame-cache-one", required_argument, 0, 'r'},
      {"frame-cache-two", required_argument, 0, 's'},
      {"calibration-cache-one", required_argument, 0, 't'},
      {"calibration-cache-two", required_argument, 0, 'u'},
      {"calibration-offset", required_argument, 0, 'v'},
      {"mask-bin", required_argument, 0, 'w'},
      {"hann-duration", required_argument, 0, 'x'},
      {"hpf-frequency", required_argument, 0, 'y'},
      {"hpf-attenuation", required_argument, 0, 'z'},
      {"hpf-order", required_argument, 0, 'A'},
      {"geo-hpf-frequency", required_argument, 0, 'B'},
      {"geo-hpf-attenuation", required_argument, 0, 'C'},
      {"geo-hpf-order", required_argument, 0, 'D'},
      {"alpha", required_argument, 0, 'E'},
      {"f-ref", required_argument, 0, 'F'},
      {"omega0", required_argument, 0, 'G'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, \
        "abd:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:" \
        "A:B:C:D:E:F:G:", long_options, &option_index);

    if (c == -1)
    {
      /* end of options, break loop */
      break;
    }

    switch(c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
        {
          break;
        }
        else
        {
          fprintf(stderr, "error parseing option %s with argument %s\n", \
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'a':
        /* help */
        display_usage();
        exit(0);
        break;

      case 'b':
        /* display version info and exit */
        fprintf(stdout, "Standalone SGWB Search Engine\n");
        XLALOutputVersionString(stderr,0);
        exit(0);
        break;

      case 'd':
        /* user tag */
        optarg_len = strlen(optarg) + 1;
        userTag = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(userTag, optarg, optarg_len);

        /* add to process_params table */
        this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                          calloc(1, sizeof(ProcessParamsTable));
        snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", prog_name);
        snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--user-tag");
        snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
        snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, "%s", optarg);
        break;

      case 'e':
        /* xml comment */
        if (strlen(optarg) > LIGOMETA_COMMENT_MAX - 1)
        {
          fprintf(stderr, "invalid argument to --%s:\n" \
              "comment must be less than %d characters\n", \
              long_options[option_index].name, LIGOMETA_COMMENT_MAX);
          exit(1);
        }
        else
        {
          snprintf(comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'f':
        /* directory for output files */
        optarg_len = strlen(optarg) + 1;
        outputPath = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputPath, optarg, optarg_len);
        if (stat(outputPath, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Directory does not exist: (%s specified)\n", \
              long_options[option_index].name, outputPath);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", outputPath);
        break;

      case 'g':
        /* start time */
        startTime = atoi(optarg);
        if (startTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        if (startTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%"LAL_INT4_FORMAT, startTime);
        break;

      case 'h':
        /* end time */
        endTime = atoi(optarg);
        if (endTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        if (endTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%"LAL_INT4_FORMAT, endTime);
        break;

      case 'i':
        /* interval duration */
        intervalDuration = atoi(optarg);
        if (intervalDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Interval duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, intervalDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", intervalDuration);
        break;

      case 'j':
        /* segment duration */
        segmentDuration = atoi(optarg);
        if (segmentDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Segment duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, segmentDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", segmentDuration);
        break;

      case 'k':
        /* resample rate */
        resampleRate = atoi(optarg);
        if (resampleRate < 2 || resampleRate > 16384 || resampleRate % 2)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Resample rate must be a power of 2 between 2 and 16384: " \
              "inclusive: (%d specified)\n", long_options[option_index].name, \
              resampleRate);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", resampleRate);

        break;

      case 'l':
        /* minimal frequency */
        fMin = atof(optarg);
        if (fMin < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Minimum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMin);
          exit(1);
        }
        /* check that min frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMin != round(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMin = round(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMin has been rounded to %f\n", fMin);
        }
        ADD_PROCESS_PARAM("float", "%e", fMin);
        break;

      case 'm':
        /* maximal frequency */
        fMax = atof(optarg);
        if (fMax < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Maximum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMax);
          exit(1);
        }
        /* check that the max frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMax != round(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMax = round(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMax has been rounded to %f\n", fMax);
        }
        ADD_PROCESS_PARAM("float", "%e", fMax);
        break;

      case 'n':
        /* ifo for first stream */
        optarg_len = strlen(optarg) + 1;
        ifoOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoOne, optarg, optarg_len);

        /* set site id for ifo one */
        if (strncmp(ifoOne, "H1", 2) == 0)
          siteOne = 0;
        else if (strncmp(ifoOne, "H2", 2) == 0)
          siteOne = 0;
        else if (strncmp(ifoOne, "L1", 2) == 0)
          siteOne = 1;
        else if (strncmp(ifoOne, "G1", 2) == 0)
          siteOne = 3;
        else
        {
          fprintf(stderr, "First IFO not recognised...\n");
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", ifoOne);
        break;

      case 'o':
        /* ifo for second stream */
        optarg_len = strlen(optarg) + 1;
        ifoTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoTwo, optarg, optarg_len);

        /* set site id for ifo two */
        if (strncmp(ifoTwo, "H1", 2) == 0)
          siteTwo = 0;
        else if (strncmp(ifoTwo, "H2", 2) == 0)
          siteTwo = 0;
        else if (strncmp(ifoTwo, "L1", 2) == 0)
          siteTwo = 1;
        else if (strncmp(ifoTwo, "G1", 2) == 0)
          siteOne = 3;
        else
        {
          fprintf(stderr, "Second IFO not recognised...\n");
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", ifoTwo);
        break;

      case 'p':
        /* channel one */
        optarg_len = strlen(optarg) + 4;
        channelOneTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(channelOneTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelOneTemp);
        break;

      case 'q':
        /* channel two */
        optarg_len = strlen(optarg) + 4;
        channelTwoTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(channelTwoTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelTwoTemp);
        break;

      case 'r':
        /* frame cache one */
        optarg_len = strlen(optarg) + 1;
        frameCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheOne, optarg, optarg_len);
        if (stat(frameCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheOne);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", frameCacheOne);
        break;

      case 's':
        /* frame cache two */
        optarg_len = strlen(optarg) + 1;
        frameCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheTwo, optarg, optarg_len);
        if (stat(frameCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheTwo);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", frameCacheTwo);
        break;

      case 't':
        /* calibration cache one */
        optarg_len = strlen(optarg) + 1;
        calCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheOne, optarg, optarg_len);
        if (stat(calCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheOne);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", calCacheOne);
        break;

      case 'u':
        /* calibration cache two */
        optarg_len = strlen(optarg) + 1;
        calCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheTwo, optarg, optarg_len);
        if (stat(calCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheTwo);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", calCacheTwo);
        break;

      case 'v':
        /* calibration time offset */
        calibOffset = atoi(optarg);
        ADD_PROCESS_PARAM("int", "%d", calibOffset);
        break;

      case 'w':
        /* number of bins to mask for frequency mask */
        maskBin = atoi(optarg);
        if (maskBin <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of bins to mask must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, maskBin);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", maskBin);
        break;

      case 'x':
        /* hann window duration */
        hannDuration = atoi(optarg);
        if (hannDuration < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Hann duartion is less than 0: (%d specified)\n", \
              long_options[option_index].name, hannDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", hannDuration);
        break;

      case 'y':
        /* high pass knee filter frequency  */
        highPassFreq = atof(optarg);
        if (highPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "High pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              highPassFreq);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassFreq);
        break;

      case 'z':
        /* high pass filter attenuation  */
        highPassAtten = atof(optarg);
        if ((highPassAtten < 0.0) || (highPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              highPassAtten);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassAtten);
        break;

      case 'A':
        /* high pass filter order  */
        highPassOrder = atoi(optarg);
        if (highPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              highPassOrder);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", highPassOrder);
        break;

      case 'B':
        /* GEO high pass knee filter frequency */
        geoHighPassFreq = atof(optarg);
        if (geoHighPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "GEO high pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassFreq);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", geoHighPassFreq);
        break;

      case 'C':
        /* GEO high pass filter attenuation */
        geoHighPassAtten = atof(optarg);
        if ((geoHighPassAtten < 0.0) || (geoHighPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassAtten);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", geoHighPassAtten);
        break;

      case 'D':
        /* GEO high pass filter order */
        geoHighPassOrder = atoi(optarg);
        if (geoHighPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              geoHighPassOrder);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", geoHighPassOrder);
        break;

      case 'E':
        /* filter spectrum exponent */
        alpha = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", alpha);
        break;

      case 'F':
        /* filter reference frequency */
        fRef = atof(optarg);
        if (fRef < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference frequency must be greater than 0: " \
              "(%f specified)\n", long_options[option_index].name, fRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", fRef);
        break;

      case 'G':
        /* filter reference omega */
        omegaRef = atof(optarg);
        if (omegaRef <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference omega_0 must be positive: (%f specified)\n", \
              long_options[option_index].name, omegaRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", omegaRef);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing options\n");
        exit(1);
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
    {
      fprintf(stderr, "%s\n", argv[optind++]);
    }
    exit(1);
  }

  /* check for required arguments */

  /* start/end time */
  if (startTime == 0)
  {
    fprintf(stderr, "--gps-start-time must be specified\n");
    exit(1);
  }
  if (endTime == 0)
  {
    fprintf(stderr, "--gps-end-time must be specified\n");
    exit(1);
  }

  /* interval duration */
  if (intervalDuration == -1)
  {
    fprintf(stderr, "--interval-duration must be specified\n");
    exit(1);
  }

  /* segment duration */
  if (segmentDuration == -1)
  {
    fprintf(stderr, "--segment-duration must be specified\n");
    exit(1);
  }

  /* min/max frequency */
  if (fMin == -1)
  {
    fprintf(stderr, "--f-min must be specified\n");
    exit(1);
  }
  if (fMax == -1)
  {
    fprintf(stderr, "--f-max must be specified\n");
    exit(1);
  }

  /* ifos */
  if (ifoOne == NULL)
  {
    fprintf(stderr, "--ifo-one must be specified\n");
    exit(1);
  }
  if (ifoTwo == NULL)
  {
    fprintf(stderr, "--ifo-two must be specified\n");
    exit(1);
  }

  /* channels */
  if (channelOne == NULL)
  {
    fprintf(stderr, "--channel-one must be specified\n");
    exit(1);
  }
  if (channelTwo == NULL)
  {
    fprintf(stderr, "--channel-two must be specified\n");
    exit(1);
  }

  /* frame cache */
  if (frameCacheOne == NULL)
  {
    fprintf(stderr, "--frame-cache-one must be specified\n");
    exit(1);
  }
  if (siteOne != siteTwo)
  {
    /* only need second frame cache if ifos differ */
    if (frameCacheTwo == NULL)
    {
      fprintf(stderr, "--frame-cache-two must be specified\n");
      exit(1);
    }
  }
  else
  {
    /* if site ids are the same then the frames for the different
     * detectors are in the same frame cache */
    frameCacheTwo = frameCacheOne;
  }

  /* calibration cache */
  if (strncmp(ifoOne, "G1", 2) != 0)
  {
    if (calCacheOne == NULL)
    {
      fprintf(stderr, "--calibration-cache-one must be specified\n");
      exit(1);
    }
  }
  if (strncmp(ifoTwo, "G1", 2) != 0)
  {
    if (calCacheTwo == NULL)
    {
      fprintf(stderr, "--calibration-cache-two must be specified\n");
      exit(1);
    }
  }

  /* calibration offset */
  if (calibOffset == -1)
  {
    fprintf(stderr, "--calibration-offset must be specified\n");
    exit(1);
  }

  /* mask */
  if (apply_mask_flag)
  {
    if (maskBin == -1)
    {
      fprintf(stderr, "--mask-bin must be specified\n");
      exit(1);
    }
  }

  /* hann duration */
  if (overlap_hann_flag)
  {
    if (hannDuration != -1)
    {
      fprintf(stderr, "Overlapping Hann windows specified, --hann-duration " \
          "will be ignored\n");
    }
  }
  else
  {
    if (hannDuration == -1)
    {
      fprintf(stderr, "--hann-duration must be specified\n");
      exit(1);
    }
  }

  /* high pass filter */
  if (high_pass_flag)
  {
    if (highPassFreq == -1)
    {
      fprintf(stderr, "--hpf-frequency must be specified\n");
      exit(1);
    }
    if (highPassAtten == -1)
    {
      fprintf(stderr, "--hpf-attenuation must be specified\n");
      exit(1);
    }
    if (highPassOrder == -1)
    {
      fprintf(stderr, "--hpf-order must be specified\n");
      exit(1);
    }
  }

  /* GEO high pass filter */
  if ((strncmp(ifoOne, "G1", 2) == 0) || (strncmp(ifoTwo, "G1", 2) == 0))
  {
    if (geoHighPassFreq == -1)
    {
      fprintf(stderr, "--geo-hpf-frequency must be specified\n");
      exit(1);
    }
    if (geoHighPassAtten == -1)
    {
      fprintf(stderr, "--geo-hpf-attenuation must be specified\n");
      exit(1);
    }
    if (geoHighPassOrder == -1)
    {
      fprintf(stderr, "--geo-hpf-order must be specified\n");
      exit(1);
    }
  }

  /* check for sensible arguments */

  /* start time same as stop time */
  if (startTime == endTime)
  {
    fprintf(stderr, "Start time same as end time; no analysis to perform\n");
    exit(1);
  }

  /* stop time before start time */
  if (startTime > endTime)
  {
    fprintf(stderr, "Invalid start/end time; end time (%d) is before " \
        "start time (%d)\n", endTime, startTime);
    exit(1);
  }

  /* interval duration must be a least 3 times the segment duration */
  if ((intervalDuration / segmentDuration) < 3)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be a least 3 times " \
        "the segment\nduration (%d)\n", intervalDuration, segmentDuration);
    exit(1);
  }

  /* interval duration must be an odd mutliple of segment duration */
  if (((intervalDuration / segmentDuration) % 2) != 1)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be an odd " \
        "multiple of the segment\nduration (%d)\n", intervalDuration, \
        segmentDuration);
    exit(1);
  }

  /* min frequency same as max */
  if (fMin == fMax)
  {
    fprintf(stderr, "Minimum frequency same as maximum; no analysis to " \
        "perform\n");
    exit(1);
  }

  /* max frequency less than min */
  if (fMin > fMax)
  {
    fprintf(stderr, "Invalid frequency band; maximum frequency (%f Hz) is " \
        "before minimum\nfrequency (%f Hz)\n", fMax, fMin);
    exit(1);
  }

  /* filter reference frequency less than min */
  if (fRef < fMin)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is less than minimum " \
        "frequency (%f Hz)\n", fRef, fMin);
    exit(1);
  }

  /* filter reference frequency greater than max */
  if (fRef > fMax)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is greater than maximum " \
        "frequency (%f Hz)\n", fRef, fMax);
    exit(1);
  }

  /* set channels */
  strcpy(channelOne, ifoOne);
  strcpy(channelTwo, ifoTwo);
  strcat(channelOne, ":");
  strcat(channelTwo, ":");
  strcat(channelOne, channelOneTemp);
  strcat(channelTwo, channelTwoTemp);
  free(channelOneTemp);
  free(channelTwoTemp);

  return;
}

/* program entry point */
INT4 main(INT4 argc, CHAR *argv[])
{
  /* lal initialisation variables */
  LALStatus status = blank_status;

  /* xml */
  CHAR baseName[FILENAME_MAX];
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;

  /* counter */
  INT4 i;

  /* data structures */
  REAL4TimeSeries *seriesOne;
  REAL4TimeSeries *seriesTwo;
  REAL4FrequencySeries *overlap;
  REAL4FrequencySeries *mask;
  REAL4FrequencySeries *omegaGW;
  REAL4Window *dataWindow;
  LIGOTimeGPS gpsStartTime;
  LIGOTimeGPS gpsEndTime;

  /* variables */
  INT4 numSegments;
  INT4 duration, durationEff, extrasec;
  INT4 segsInInt;
  INT4 segmentLength;
  INT4 padData;
  REAL8 deltaF;
  INT4 filterLength;
  INT4 numFMin;
  INT4 numFMax;

  /* error handler */
  status.statusPtr = NULL;
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc(1, sizeof(ProcessTable));
  XLALGPSTimeNow(&proctable.processTable->start_time);
  XLALPopulateProcessTable(proctable.processTable, prog_name, \
      LALAPPS_VCS_IDENT_ID, LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) \
                    calloc(1, sizeof(ProcessParamsTable));
  memset(comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR));

  /* parse command line options */
  parse_options(argc, argv);

  /* get xml file basename */
  if (userTag)
  {
    snprintf(baseName, FILENAME_MAX, "%s%s-STOCHASTIC_%s_%d-%d", \
        ifoOne, ifoTwo, userTag, startTime, (endTime - startTime));
  }
  else
  {
    snprintf(baseName, FILENAME_MAX, "%s%s-STOCHASTIC-%d-%d", \
        ifoOne, ifoTwo, startTime, (endTime - startTime));
  }

  /* get number of segments */
  duration = endTime - startTime;
  numSegments = duration / segmentDuration;
  segsInInt = intervalDuration / segmentDuration;

  /* recentre */
  if (recentre_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Recentring within data stream...\n");

    durationEff = numSegments * segmentDuration;
    extrasec = duration - durationEff;
    startTime += extrasec / 2;
    endTime = startTime + durationEff;
  }

  /* add a resample buffer, if required */
  if ((resampleRate) || (high_pass_flag))
    padData = 1;
  else
    padData = 0;

  /* initialise gps time structures */
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0;
  gpsEndTime.gpsSeconds = endTime;
  gpsEndTime.gpsNanoSeconds = 0;

  /* read data */
  seriesOne = get_time_series(ifoOne, frameCacheOne, channelOne, \
      gpsStartTime, gpsEndTime, resampleRate, highPassOrder, highPassFreq, \
      highPassAtten, geoHighPassOrder, geoHighPassFreq, geoHighPassAtten, \
      padData);
  seriesTwo = get_time_series(ifoTwo, frameCacheTwo, channelTwo, \
      gpsStartTime, gpsEndTime, resampleRate, highPassOrder, highPassFreq, \
      highPassAtten, geoHighPassOrder, geoHighPassFreq, geoHighPassAtten, \
      padData);

  /* check that the two series have the same sample rate */
  if (seriesOne->deltaT != seriesTwo->deltaT)
  {
    fprintf(stderr, "Series have different sample rates...\n");
    exit(1);
  }
  else
  {
    /* get resample rate, if required */
    if (!resampleRate)
      resampleRate = (INT4)(1./seriesOne->deltaT);
  }

  /* get deltaF for optimal filter */
  deltaF = 1./(REAL8)PSD_WINDOW_DURATION;

  /* get bins for min and max frequencies */
  numFMin = (INT4)(fMin / deltaF);
  numFMax = (INT4)(fMax / deltaF);

  /* get lengths */
  filterLength = numFMax - numFMin + 1;
  segmentLength = segmentDuration * resampleRate;

  if (vrbflg)
    fprintf(stdout, "Generating data segment window...\n");

  /* for overlapping hann windows, the hann window length is the segment
   * length */
  if (overlap_hann_flag)
    hannDuration = segmentDuration;

  /* create window for data */
  dataWindow = data_window(seriesOne->deltaT, segmentLength, hannDuration);

  /* generate overlap reduction function */
  if (vrbflg)
    fprintf(stdout, "Generating the overlap reduction function...\n");
  overlap = overlap_reduction_function(&status, filterLength, fMin, deltaF, \
      siteOne, siteTwo);

  /* generate omegaGW */
  if (vrbflg)
    fprintf(stdout, "Generating spectrum for optimal filter...\n");
  omegaGW = omega_gw(&status, alpha, fRef, omegaRef, filterLength, \
      fMin, deltaF);

  /* frequency mask */
  if (apply_mask_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Applying frequency mask to spectrum..\n");

    /* generate frequency mask */
    mask = frequency_mask(fMin, deltaF, filterLength, maskBin);

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
      omegaGW->data->data[i] *= mask->data->data[i];

    /* destroy frequency mask */
    XLALDestroyREAL4FrequencySeries(mask);
  }

  /* perform search */
  stochHead = stochastic_search(&status, seriesOne, seriesTwo, overlap, \
      omegaGW, dataWindow, numSegments, filterLength, segmentDuration, \
      segsInInt, segmentLength, PSD_WINDOW_DURATION, ifoOne, ifoTwo, \
      channelOne, channelTwo, calCacheOne, calCacheTwo, calibOffset, \
      fMin, fMax, fRef);

  /* save out xml table */
  save_xml_file(&status, prog_name, outputPath, baseName, \
      stochHead, proctable, procparams, this_proc_param, comment);

  /* cleanup */
  XLALDestroyREAL4FrequencySeries(overlap);
  XLALDestroyREAL4FrequencySeries(omegaGW);
  XLALDestroyREAL4Window(dataWindow);
  XLALDestroyREAL4TimeSeries(seriesOne);
  XLALDestroyREAL4TimeSeries(seriesTwo);

  /* free memory used in the stochastic xml table */
  while(stochHead)
  {
    thisStoch = stochHead;
    stochHead = stochHead->next;
    LALFree(thisStoch);
  }

  /* free calloc'd memory */
  if (strcmp(frameCacheOne, frameCacheTwo))
  {
    free(frameCacheOne);
    free(frameCacheTwo);
  }
  else
  {
    free(frameCacheOne);
  }
  free(calCacheOne);
  free(calCacheTwo);
  free(channelOne);
  free(channelTwo);
  free(ifoOne);
  free(ifoTwo);
  free(userTag);
  free(outputPath);

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit(0);
}

/*
 * vim: et
 */
