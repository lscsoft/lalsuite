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
#include "skygrid.h"
#include <string.h>

int lalDebugLevel = 0;
int verbosity_level = 0;
struct gengetopt_args_info args_info;


int
main(int argc, char **argv)
{
  static LALStatus          status;
  EphemerisData             ephemeris_data;
  skygrid_t grid_cros_sq = NULL;
  skygrid_t grid_plus_sq = NULL;
  skygrid_t grid_sum_sq  = NULL;
  skygrid_t grid_relfreq = NULL;
  
  status.statusPtr = NULL;
  
  /* parse command line options */
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.verbosity_given)
    verbosity_level = args_info.verbosity_arg;

  if (args_info.debug_given)
    lalDebugLevel = args_info.debug_arg;

  if (!args_info.count_ra_given)
  	args_info.count_ra_arg = args_info.n_ra_arg-args_info.start_ra_arg;
    
  if (!args_info.count_dec_given)
  	args_info.count_dec_arg = args_info.n_dec_arg-args_info.start_dec_arg;
    
  if(lalDebugLevel!=0)printf("lalDebugLevel = %d\n", lalDebugLevel);
    
  /* initalize ephemeris_data */
  init_ephemeris(&status, &ephemeris_data);

  /* initialize skygrid_t stuff - must be after cmdline is parsed */
  init_skygrid(&status);

  /* allocate skygrids */
  (void *)alloc_skygrid(&status, &grid_cros_sq);
  (void *)alloc_skygrid(&status, &grid_plus_sq);
  (void *)alloc_skygrid(&status, &grid_sum_sq);
  (void *)alloc_skygrid(&status, &grid_relfreq);
  
  if (lalDebugLevel)
  {
    printf("grid_cros_sq->dimLength->data = %d x %d\n", 
           grid_cros_sq->dimLength->data[0], grid_cros_sq->dimLength->data[1]);
    printf("grid_plus_sq->dimLength->data = %d x %d\n", 
           grid_plus_sq->dimLength->data[0], grid_plus_sq->dimLength->data[1]);
    printf("grid_sum_sq->dimLength->data = %d x %d\n", 
           grid_sum_sq->dimLength->data[0], grid_sum_sq->dimLength->data[1]);
    printf("grid_relfreq->dimLength->data = %d x %d\n", 
           grid_relfreq->dimLength->data[0], grid_relfreq->dimLength->data[1]);
  }

  if (args_info.single_source_given)
    generate_timeseries_response(&status);
  else if (args_info.whole_sky_given)
    compute_skygrid(&status, &ephemeris_data, grid_cros_sq,
                    grid_plus_sq, grid_sum_sq, grid_relfreq);
    
  
  /*
   * Housekeeping
   */
  free_skygrid(&status, &grid_cros_sq);
  free_skygrid(&status, &grid_plus_sq);
  free_skygrid(&status, &grid_sum_sq);
  free_skygrid(&status, &grid_relfreq);
  
  cleanup_skygrid(&status);
    
  cleanup_ephemeris(&status, &ephemeris_data);
  
  LALCheckMemoryLeaks();

  return 0;
}

