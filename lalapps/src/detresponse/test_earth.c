/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

#include "config.h"
#include <lalapps.h>

#include "cmdline.h"
#include "skygrid.h"
#include "util.h"
#include <string.h>

int lalDebugLevel = 7;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int
main(int argc, char **argv)
{
  static LALStatus s;
  EphemerisData ephemerides;
  EarthState    earth;
  REAL8         earth_phi; /* azi. position of Earth in solar
                              system barycenter */
  LIGOTimeGPS   gps;
  UINT4         i;

  s.statusPtr = NULL;
    
  init_ephemeris(&s, &ephemerides);
  
  gps.gpsSeconds = 669130273; /* vernal equinox 2001 */
  
  for (i = 0; i < 185; ++i)
  {
    gps.gpsSeconds += i * 3600*24;
    LALBarycenterEarth(&s, &earth, &gps, &ephemerides);
    printf("%d {% e, % e, % e} {% e, % e, % e}\n",
           gps.gpsSeconds,
           earth.posNow[0], earth.posNow[1], earth.posNow[2],
           earth.velNow[0], earth.velNow[1], earth.velNow[2]);  
  }
  
  cleanup_ephemeris(&s, &ephemerides);
  
  return 0;
}