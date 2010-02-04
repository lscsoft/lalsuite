/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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

extern int lalDebugLevel;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int
main(int argc, char **argv)
{
  static LALStatus          status;
  EphemerisData             *ephemeris_data;
  skygrid_t grid_cros_sq = NULL;
  skygrid_t grid_plus_sq = NULL;
  skygrid_t grid_sum_sq  = NULL;
  skygrid_t grid_relfreq = NULL;

  lalDebugLevel = 0;

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
  if ( (ephemeris_data = InitEphemeris (NULL, "03-06")) == NULL ) {
    XLALPrintError ("Failed to initialized '03-06' ephemeris, make sure to set LAL_DATA_PATH correctly!\n");
    exit(1);
  }

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
    compute_skygrid(&status, ephemeris_data, grid_cros_sq,
                    grid_plus_sq, grid_sum_sq, grid_relfreq);
    
  
  /*
   * Housekeeping
   */
  free_skygrid(&status, &grid_cros_sq);
  free_skygrid(&status, &grid_plus_sq);
  free_skygrid(&status, &grid_sum_sq);
  free_skygrid(&status, &grid_relfreq);
  
  cleanup_skygrid(&status);

  XLALFree ( ephemeris_data->ephemE );
  XLALFree ( ephemeris_data->ephemS );
  XLALFree ( ephemeris_data );

  
  LALCheckMemoryLeaks();

  return 0;
}

