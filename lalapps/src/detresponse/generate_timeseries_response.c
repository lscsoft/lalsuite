/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include "config.h"
#include "detresponse.h"

/*
 * TODO: make file names depend on name of source, and of detector
 */

void generate_timeseries_response(LALStatus * status)
{
  LALSource                 source;
  LALFrDetector             frdetector;
  LALDetector               detector;
  LALDetAndSource           det_and_src = {NULL,NULL};
  LALDetAMResponseSeries    det_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  REAL4TimeSeries           tmp_series;
  LALTimeIntervalAndNSample time_info;
  CHAR  cross_file_name[LALNameLength];
  CHAR  plus_file_name[LALNameLength];
  CHAR  sum_file_name[LALNameLength];

  

  /* null out strings */
  cross_file_name[0] = '\0';
  plus_file_name[0] = '\0';
  sum_file_name[0] = '\0';

  if (args_info.detector_given)
    {
      if (mystrncasecmp(args_info.detector_arg,"lho", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"llo", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"virgo",LALNameLength)== 0)
        detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"geo", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"tama",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"cit", LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"test", LALNameLength) == 0)
        {
          set_detector_params(status, &frdetector, &detector,
                              "TEST - North Pole",
                              0., LAL_PI_2, 0.,
                              0., LAL_PI_2,
                              0., 0.);
        }
      else
        {
          /* unknown detector -- exit with error */
          fprintf(stderr, "lalapps_detresponse: ERROR: Unknown detector '%s'\n",
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

  (void)mystrlcpy(cross_file_name, args_info.source_name_arg, LALNameLength);
  (void)strncat(cross_file_name, "_cross.txt",
                (LALNameLength - strlen(args_info.source_name_arg)));
  (void)mystrlcpy(plus_file_name, args_info.source_name_arg, LALNameLength);
  (void)strncat(plus_file_name, "_plus.txt",
                (LALNameLength - strlen(args_info.source_name_arg)));
  (void)mystrlcpy(sum_file_name, args_info.source_name_arg, LALNameLength);
  (void)strncat(sum_file_name, "_sum.txt",
                (LALNameLength - strlen(args_info.source_name_arg)));
    
  if ((lalDebugLevel & 1) && (verbosity_level & 4))
    {
      printf("cross_file_name = %s\n", cross_file_name);
      printf("plus_file_name  = %s\n", plus_file_name);
      printf("sum_file_name   = %s\n", sum_file_name);
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
  else
    {
      fprintf(stderr, "ERROR: time info incomplete\n");
      exit(5);
    }

  if (verbosity_level & 4)
    print_time_info(&time_info);


  plus_series.data   = NULL;
  cross_series.data  = NULL;
  scalar_series.data = NULL;
  tmp_series.data    = NULL;

  det_response_series.pPlus   = &(plus_series);
  det_response_series.pCross  = &(cross_series);
  det_response_series.pScalar = &(scalar_series);

  LALSCreateVector(status, &(det_response_series.pPlus->data),
                   time_info.nSample);
  LALSCreateVector(status, &(det_response_series.pCross->data),
                   time_info.nSample);
  LALSCreateVector(status, &(det_response_series.pScalar->data),
                   time_info.nSample);
  LALSCreateVector(status, &(tmp_series.data), time_info.nSample);

  LALComputeDetAMResponseSeries(status, &det_response_series,
                                &det_and_src, &time_info);


  /* want F+^2, Fx^2, and (F+^2 + Fx^2) */
  square_timeseries(det_response_series.pPlus);
  square_timeseries(det_response_series.pCross);
  add_timeseries(&tmp_series, det_response_series.pPlus,
                 det_response_series.pCross);

  LALSPrintTimeSeries(det_response_series.pPlus, plus_file_name);
  LALSPrintTimeSeries(det_response_series.pCross, cross_file_name);
  LALSPrintTimeSeries(&tmp_series, sum_file_name);

  /*
   * house keeping
   */
  LALSDestroyVector(status, &(det_response_series.pPlus->data));
  LALSDestroyVector(status, &(det_response_series.pCross->data));
  LALSDestroyVector(status, &(det_response_series.pScalar->data));
  LALSDestroyVector(status, &(tmp_series.data));

  LALCheckMemoryLeaks();
}
