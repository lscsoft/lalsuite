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
  UINT4            num_ra, num_dec;
  gridding_t       g;
  LIGOTimeGPS      gps;
  LALLeapSecAccuracy acc = LALLEAPSEC_LOOSE;

  s.statusPtr = NULL;
  
  LALGPSTimeNow (&s, &gps, &acc); 
  
  printf("RUN 1\n");
  
  init_gridding(&g);
  
  num_ra = 24;
  num_dec = 11;
  make_gridding(&s, &g, num_ra, DETRESP_REGGRID, num_dec, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, (char *)NULL);
  
  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 2 - Autumnal Equi. 2003\n");
  
  gps.gpsSeconds = 748349233; /* Autumnal Equi. 2003 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, num_dec, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, "autumn2003.dat");

  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 3 - Winter Sols. 2003\n");
  
  gps.gpsSeconds = 756111853; /* Winter Sols. 2003 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, num_dec, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, "winter2003.dat");

  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 4 - Vernal Equi. 2004\n");
  
  gps.gpsSeconds = 763800553; /* Vernal Equi. 2004 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, num_dec, DETRESP_REGGRID, &gps);
  
  print_gridding(&g, "spring2004.dat");

  cleanup_gridding(&s, &g);

  LALCheckMemoryLeaks();

  return 0;
}
