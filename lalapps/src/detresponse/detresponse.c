/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

/*
 * Computes detresponse for a given source and given detector, for a given
 * time period.
 *
 * Inputs: source params - RA, dec, orientation
 *         IFO params - name, mode
 *         time info - start time, duration, interval
 * 
 */

#include "config.h"
#include <lalapps.h>

#include "detresponse.h"
#include <string.h>

#include "cmdline.h"
#include "strncasecmp.h"

static int lalDebugLevel = 0;
static int verbosity_level = 0;

int
main(int argc, char **argv)
{
  struct gengetopt_args_info args_info;
  
  static LALStatus          status;
  LALSource                 source;
  LALDetector               detector;
  LALDetAndSource           det_and_src = {NULL,NULL};
  LALDetAMResponseSeries    det_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  LALTimeIntervalAndNSample time_info;
  char  cross_file_name[LALNameLength];
  char  plus_file_name[LALNameLength];

  /* null out strings */
  cross_file_name[0] = '\0';
  plus_file_name[0] = '\0';
  
  /* cmdline parser */
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.verbosity_given)
    verbosity_level = args_info.verbosity_arg;

  if (args_info.debug_given)
    lalDebugLevel = args_info.debug_arg;

  if (args_info.detector_given)
    {
      if (strncasecmp(args_info.detector_arg, "lho", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
      else if (strncasecmp(args_info.detector_arg, "llo", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
      else if (strncasecmp(args_info.detector_arg, "virgo", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      else if (strncasecmp(args_info.detector_arg, "geo", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      else if (strncasecmp(args_info.detector_arg, "tama", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      else if (strncasecmp(args_info.detector_arg, "cit", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
      else
        {
          /* unknown detector -- exit with error */
          fprintf(stderr, "lalapps_detresponse: Unknown detector '%s'\n",
                  args_info.detector_arg);
          exit(2);
        }
    }

  if (verbosity_level & 4)
    PrintLALDetector(&detector);

  if (args_info.right_ascension_given &&
      args_info.declination_given &&
      args_info.orientation_given)
    {
      set_source_params(&source, args_info.source_name_arg,
                        args_info.right_ascension_arg,
                        args_info.declination_arg,
                        args_info.orientation_arg);
    }

  if (verbosity_level & 4)
    print_source(&source);

  (void)strlcpy(cross_file_name, args_info.source_name_arg, LALNameLength);
  (void)strncat(cross_file_name, "_cross.txt",
                (LALNameLength - strlen(args_info.source_name_arg)));
  (void)strlcpy(plus_file_name, args_info.source_name_arg, LALNameLength);
  (void)strncat(plus_file_name, "_plus.txt",
                (LALNameLength - strlen(args_info.source_name_arg)));
  
  if ((lalDebugLevel & 1) && (verbosity_level & 4))
    {
      printf("cross_file_name = %s\n", cross_file_name);
      printf("plus_file_name  = %s\n", plus_file_name);
    }

  det_and_src.pDetector = &detector;
  det_and_src.pSource   = &source;


  if (args_info.start_time_sec_given &&
      args_info.nsample_given &&
      args_info.sampling_interval_given)
    {
      time_info.accuracy             = LALLEAPSEC_STRICT;
      
      if (verbosity_level == 13)
        {
          time_info.epoch.gpsSeconds     = 709398013;
          time_info.epoch.gpsNanoSeconds =         0;
          time_info.deltaT               =       864;
          time_info.nSample              =       100;
        }
      else
        {
          time_info.epoch.gpsSeconds     = args_info.start_time_sec_arg;
          time_info.epoch.gpsNanoSeconds = args_info.start_time_nanosec_arg;
          time_info.deltaT               = args_info.sampling_interval_arg;
          time_info.nSample              = args_info.nsample_arg;
        }
    }

  if (verbosity_level & 4)
    print_time_info(&time_info);


  plus_series.data   = NULL;
  cross_series.data  = NULL;
  scalar_series.data = NULL;

  det_response_series.pPlus   = &(plus_series);
  det_response_series.pCross  = &(cross_series);
  det_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(det_response_series.pPlus->data), 1);
  LALSCreateVector(&status, &(det_response_series.pCross->data), 1);
  LALSCreateVector(&status, &(det_response_series.pScalar->data), 1);

  LALComputeDetAMResponseSeries(&status, &det_response_series,
                                          &det_and_src, &time_info);

  LALSPrintTimeSeries(det_response_series.pPlus, plus_file_name);
  LALSPrintTimeSeries(det_response_series.pCross, cross_file_name);


  /*
   * house keeping
   */
  LALSDestroyVector(&status, &(det_response_series.pPlus->data));
  LALSDestroyVector(&status, &(det_response_series.pCross->data));
  LALSDestroyVector(&status, &(det_response_series.pScalar->data));

  return 0;
}
