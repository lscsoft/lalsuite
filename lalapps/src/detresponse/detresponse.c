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

int lalDebugLevel = 0;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int
main(int argc, char **argv)
{
  static LALStatus          status;

  /* parse command line options */
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.verbosity_given)
    verbosity_level = args_info.verbosity_arg;

  if (args_info.debug_given)
    lalDebugLevel = args_info.debug_arg;

  generate_timeseries_response(&status);

  return 0;
}

