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
#include <lalapps.h>

#include "cmdline.h"
#include "make_gridding.h"
#include "skygrid.h"
#include "util.h"
#include <string.h>

int lalDebugLevel = 7;
int verbosity_level = 0;
struct gengetopt_args_info args_info;

int
main(int argc, char **argv)
{
  static LALStatus   s;
  UINT4              num_ra, num_dec;
  gridding_t         g;
  LIGOTimeGPS        gps;
  LALLeapSecAccuracy acc = LALLEAPSEC_LOOSE;
  EphemerisData      ephem;
  INT4               leap_secs;

  s.statusPtr = NULL;
  
  XLALGPSTimeNow (&gps); 
  init_ephemeris(&s, &ephem);

  
  printf("RUN 1\n");
  
  init_gridding(&g);
  
  num_ra = 24;
  num_dec = 11;
  make_gridding(&s, &g, num_ra, DETRESP_REGGRID, 
                num_dec, DETRESP_REGGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "reg_dec_reg_ra.txt", DETRESP_HUMANREAD);
  print_gridding(&g, "reg_dec_reg_ra.dat", DETRESP_XYPAIRS_ASCII);
  
  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 2 - Autumnal Equi. 2003\n");
  
  gps.gpsSeconds = 748349233; /* Autumnal Equi. 2003 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, 
                num_dec, DETRESP_REGGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "autumn2003.txt", DETRESP_HUMANREAD);
  print_gridding(&g, "autumn2003.dat", DETRESP_XYPAIRS_ASCII);

  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 3 - Winter Sols. 2003\n");
  
  gps.gpsSeconds = 756111853; /* Winter Sols. 2003 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, 
                num_dec, DETRESP_REGGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "winter2003.txt", DETRESP_HUMANREAD);
  print_gridding(&g, "winter2003.dat", DETRESP_XYPAIRS_ASCII);

  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");
  
  printf("RUN 4 - Vernal Equi. 2004\n");
  
  gps.gpsSeconds = 763800553; /* Vernal Equi. 2004 */
  
  num_ra = 100;
  num_dec = 51;
  make_gridding(&s, &g, num_ra, DETRESP_IRRGRID, 
                num_dec, DETRESP_REGGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "spring2004.txt", DETRESP_HUMANREAD);
  print_gridding(&g, "spring2004.dat", DETRESP_XYPAIRS_ASCII);

  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");

  printf("RUN 5 - regular Dec, variable RA\n");
  
  init_gridding(&g);
  
  num_ra = 48;
  num_dec = 22;
  make_gridding(&s, &g, num_ra, DETRESP_VARGRID, 
                num_dec, DETRESP_REGGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "reg_dec_var_ra.dat", DETRESP_XYPAIRS_ASCII);
  
  cleanup_gridding(&s, &g);
  
  printf("\n*   *   *   *   *   *   *   *   *   *\n");

  printf("RUN 6 - irregular Dec, variable RA\n");
  
  init_gridding(&g);
  
  num_ra = 48;
  num_dec = 22;
  make_gridding(&s, &g, num_ra, DETRESP_VARGRID, 
                num_dec, DETRESP_IRRGRID, &ephem, &gps, acc);
  
  print_gridding(&g, "irr_dec_var_ra.txt", DETRESP_HUMANREAD);
  print_gridding(&g, "irr_dec_var_ra.dat", DETRESP_XYPAIRS_ASCII);

  cleanup_gridding(&s, &g);

  printf("\n*   *   *   *   *   *   *   *   *   *\n");

  printf("RUN 7 - RA and Dec in separate files\n");

  init_gridding(&g);

  num_ra = 48;
  num_dec = 22;
  make_gridding(&s, &g, num_ra, DETRESP_REGGRID,
		num_dec, DETRESP_REGGRID, &ephem, &gps, acc);

  print_ra_grid(&g, "test_ra_reg.txt");
  print_dec_grid(&g, "test_dec_reg.txt");

  cleanup_gridding(&s, &g);
  
  /*
   * Housekeeping
   */
  cleanup_ephemeris(&s, &ephem);

  LALCheckMemoryLeaks();

  return 0;
}
