/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

#include "config.h"
#include <lalapps.h>

#include "cmdline.h"
#include "make_gridding.h"
#include <string.h>

int lalDebugLevel = 7;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int
main(int argc, char **argv)
{
  static LALStatus s;
  gridding_t       g;
  LIGOTimeGPS      gps;
  LALLeapSecAccuracy acc = LALLEAPSEC_LOOSE;

  s.statusPtr = NULL;
  
  LALGPSTimeNow (&s, &gps, &acc); 
  
  printf("RUN 1\n");
  
  init_gridding(&g);
  
  make_gridding(&s, &g, 24, DETRESP_REGGRID, 24, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, (char *)NULL);
  
  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 2\n");
  
  init_gridding(&g);
  
  make_gridding(&s, &g, 24, DETRESP_IRRGRID, 24, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, (char *)NULL);

  cleanup_gridding(&s, &g);

  LALCheckMemoryLeaks();

  return 0;
}
