/*
*  Copyright (C) 2007 David Chin
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

#include "config.h"

#include <lal/DetectorSite.h>

#include <lalapps.h>

#include "set_source_params.h"
#include "util.h"
#include "skygrid.h"
#include "make_gridding.h"

#include <string.h>

int lalDebugLevel = 0;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int main(int argc, char **argv)
{
  static LALStatus s;
  FILE            *cross = NULL;
  FILE            *plus  = NULL;
  LALGPSandAcc     time_info;
  EphemerisData    ephem;
  LALSource        source;
  LALFrDetector    frdetector;
  LALDetector      detector;
  LALDetAndSource  det_and_src = {NULL,NULL};
  LALDetAMResponse response;
  gridding_t       g;
  UINT4            i, j;
  UINT4            n_ra, n_dec;
  REAL8		   src_orientation;
  
  if (argc != 7)
  {
    fprintf(stderr, "Error: need N_RA, N_DEC, bar longitude, bar latitude, bar azimuth, src orientation\n");
    exit(13);
  }

  /* 
   * Initializations
   */
  s.statusPtr = NULL;
  init_ephemeris(&s, &ephem);
  init_gridding(&g);

  n_ra = atoi(argv[1]);
  n_dec = atoi(argv[2]);

  cross = xfopen("bar_cross.txt","wo");
  plus  = xfopen("bar_plus.txt","wo");

  src_orientation = atof(argv[6]);

  /*
   * there are no cached bar detectors, so need to construct one
   */
  (void)mystrlcpy(frdetector.name, "FOOBAR", LALNameLength);
  frdetector.vertexLongitudeRadians = atof(argv[3]);
  frdetector.vertexLatitudeRadians  = atof(argv[4]);
  frdetector.vertexElevation        = 0.0; 
  frdetector.xArmAltitudeRadians    = 0.0;
  frdetector.xArmAzimuthRadians     = atof(argv[5]);
  frdetector.yArmAltitudeRadians    = 0.0;
  frdetector.yArmAzimuthRadians     = 0.0;

  LALCreateDetector(&s, &detector, &frdetector, LALDETECTORTYPE_CYLBAR);

  if (verbosity_level > 0)
  {
    PrintLALDetector(&detector);
    printf("\n\n");
  }
  
  /*
   * set up time and sampling
   */
  time_info.accuracy           = LALLEAPSEC_STRICT;
  time_info.gps.gpsSeconds     = 709398013;
  time_info.gps.gpsNanoSeconds =         0;

  make_gridding(&s, &g, n_ra, DETRESP_REGGRID, n_dec, DETRESP_REGGRID,
                &ephem, &(time_info.gps));

  print_ra_grid(&g, "bar_ra_grid.txt");
  print_dec_grid(&g, "bar_dec_grid.txt");

  for (i = 0; i < g.ra->length; ++i)
  {
    for (j = 0; j < g.dec->length; ++j)
    {
      /* 
       * set up source
       */
      set_source_params(&source, "BIGSTAR", g.ra->data[i],
                        g.dec->data[j], src_orientation);
      
      /* print_source(&source); */
    
      /*
       * set up input structures
       */
      det_and_src.pDetector = &detector;
      det_and_src.pSource   = &source;
    
      LALComputeDetAMResponse(&s, &response, &det_and_src, &time_info);

      if (verbosity_level > 0)
      {
	print_response(&response);
	printf("\n\n");
      }

      fprintf(cross, "%14e\t", response.cross);
      fprintf(plus, "%14e\t", response.plus);
    }
    fprintf(cross, "\n");
    fprintf(plus, "\n");
  }

  /*
   * housekeeping
   */
  fclose(cross);
  fclose(plus);

  cleanup_ephemeris(&s, &ephem);
  cleanup_gridding(&s, &g);
            
  LALCheckMemoryLeaks();
  
  return 0;
}
