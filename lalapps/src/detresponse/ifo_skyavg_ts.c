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

#define DEGTORAD (LAL_PI/180.)

int lalDebugLevel = 0;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int main(int argc, char **argv)
{
  static LALStatus s;
  FILE            *cross = NULL;
  FILE            *plus  = NULL;
  FILE            *times = NULL;
  LALGPSandAcc     time_info, delta_t;
  EphemerisData    ephem;
  LALSource        source;
  LALFrDetector    frdetector;
  LALDetector      detector;
  LALDetAndSource  det_and_src = {NULL,NULL};
  LALDetAMResponse response;
  gridding_t       g;
  gridding_t       delta_g;
  UINT4            i, j, k;
  UINT4            n_ra, n_dec, n_timestep;
  REAL8		   src_orientation;
  REAL8            avg_plus, avg_cross;
  
  
  if (argc != 10)
  {
    fprintf(stderr, "Error: need N_RA, N_DEC, ifo long., ifo lat., x-arm azi, src orien, start time, sampling interval, no. of timesteps.\n");
    fprintf(stderr, "       all angles in decimal degrees.\n");
    exit(13);
  }

  /* 
   * Initializations
   */
  s.statusPtr = NULL;
  init_ephemeris(&s, &ephem);
  init_gridding(&g);
  init_gridding(&delta_g);

  n_ra = atoi(argv[1]);
  n_dec = atoi(argv[2]);

  cross = xfopen("ifo_cross_skyavg_ts.txt","wo");
  plus  = xfopen("ifo_plus_skyavg_ts.txt","wo");
  times = xfopen("times.txt","wo");

  avg_plus = 0.;
  avg_cross = 0.;

  /*
   * Create dummy IFO
   */
  (void)mystrlcpy(frdetector.name, "FOOBAR", LALNameLength);
  frdetector.vertexLongitudeRadians = atof(argv[3]) * DEGTORAD;
  frdetector.vertexLatitudeRadians  = atof(argv[4]) * DEGTORAD;
  frdetector.vertexElevation        = 0.0; 
  frdetector.xArmAltitudeRadians    = 0.0;
  frdetector.xArmAzimuthRadians     = atof(argv[5]) * DEGTORAD;
  frdetector.yArmAltitudeRadians    = 0.0;
  frdetector.yArmAzimuthRadians     = frdetector.xArmAzimuthRadians + LAL_PI/2.;

  src_orientation = atof(argv[6]) * DEGTORAD;

  LALCreateDetector(&s, &detector, &frdetector, LALDETECTORTYPE_IFODIFF);

  if (verbosity_level > 0)
  {
    PrintLALDetector(&detector);
    printf("\n\n");
  }
  
  /*
   * set up time and sampling
   */
  time_info.accuracy           = LALLEAPSEC_STRICT;
  time_info.gps.gpsSeconds     = atoi(argv[7]);  /* start time */
  time_info.gps.gpsNanoSeconds =         0;

  delta_t.accuracy           = LALLEAPSEC_STRICT;
  delta_t.gps.gpsSeconds     = atoi(argv[8]);  /* sampling interval */
  delta_t.gps.gpsNanoSeconds =         0;

  n_timestep = atoi(argv[9]); /* no. of timesteps */

  /* here's the thing: the average over the whole sky is the
   * integral of the antenna pattern divided by 4pi 
   * i.e. Int f(theta,phi) sin(theta) d_theta d_phi
   * So, need a grid of the delta-theta and -phi
   */

  /* the gridding doesn't change with time */
  make_gridding(&s, &g, n_ra, DETRESP_REGGRID, n_dec, DETRESP_REGGRID,
                &ephem, &(time_info.gps), time_info.accuracy);
  make_gridding(&s, &delta_g, n_ra, DETRESP_REGGRID, n_dec, DETRESP_REGGRID,
                &ephem, &(time_info.gps), time_info.accuracy);

  for (i = 0; i < g.ra->length - 1; ++i)
  {
    delta_g.ra->data[i] = g.ra->data[i+1] - g.ra->data[i];
  }

  for (i = 0; i < g.dec->length - 1; ++i)
  {
    delta_g.dec->data[i] = g.dec->data[i+1] - g.dec->data[i];
  }

  for (k = 0; k < n_timestep; ++k) 
  {
    avg_plus  = 0.;
    avg_cross = 0.;

    for (i = 0; i < g.ra->length - 1; ++i)
    {
      for (j = 0; j < g.dec->length - 1; ++j)
      {
        /* 
         * set up source
         */
        set_source_params(&source, "BIGSTAR", g.ra->data[i],
                          g.dec->data[j], src_orientation);
      
        /*
         * set up input structures
         */
        det_and_src.pDetector = &detector;
        det_and_src.pSource   = &source;
      
        LALComputeDetAMResponse(&s, &response, &det_and_src, &time_info);

        avg_plus  += (response.plus * response.plus) * sin(LAL_PI_2 - g.dec->data[j]) * delta_g.ra->data[i] * delta_g.dec->data[j];
        avg_cross += (response.cross * response.cross) * sin(LAL_PI_2 - g.dec->data[j]) * delta_g.ra->data[i] * delta_g.dec->data[j];
      }
    }

    printf("% .4g\t%d\n", avg_plus, g.ra->length * g.dec->length);

    avg_plus  /= LAL_TWOPI * 2.;
    avg_cross /= LAL_TWOPI * 2.;

    fprintf(plus,  "% .14e\n", avg_plus);
    fprintf(cross, "% .14e\n", avg_cross);
    fprintf(times, "%ld\n", time_info.gps.gpsSeconds);

    time_info.gps.gpsSeconds += delta_t.gps.gpsSeconds;
  }

  /*
   * housekeeping
   */
  fclose(cross);
  fclose(plus);
  fclose(times);

  cleanup_ephemeris(&s, &ephem);
  cleanup_gridding(&s, &g);
  cleanup_gridding(&s, &delta_g);
            
  LALCheckMemoryLeaks();
  
  return 0;
}
