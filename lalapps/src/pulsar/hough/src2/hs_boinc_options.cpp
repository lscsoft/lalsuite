/*
*  Copyright (C) 2007 Bernd Machenschalk
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

/* C++ -> C interface for BOINC_OPTION handling etc.
*/

/*
#ifdef _MSC_VER
extern ostream cout;
extern ostream cerr;
#else
using namespace std;
#include <iostream>
#endif
*/

#include "hs_boinc_options.h"
#include "boinc_api.h"

char* rcsid = "$Id$";
BOINC_OPTIONS eah_boinc_options;
APP_INIT_DATA eah_app_init_data;

/* BOINC APIv6 stuff */

double boincv6_skypos_rac = 0;
double boincv6_skypos_dec = 0;
double boincv6_fraction_done = 0;

#ifdef BOINC_APIV6
#include "graphics2.h"

HS_SHMEM* shmem = NULL;

static void update_shmem(void) {

  if (!shmem) return;

  boinc_get_init_data(eah_app_init_data);

  shmem->skypos_rac    = boincv6_skypos_rac;
  shmem->skypos_dec    = boincv6_skypos_dec;
  shmem->fraction_done = boinc_get_fraction_done(); // boincv6_fraction_done;
  shmem->cpu_time      = boinc_worker_thread_cpu_time();;
  boinc_get_status(&shmem->status);

  snprintf(shmem->xml, sizeof(shmem->xml),
	   "<graphics_info>\n"
	   "  <skypos_rac>%f</skypos_rac>\n"
	   "  <skypos_dec>%f</skypos_dec>\n"
	   "  <fraction_done>%f</fraction_done>\n"
	   "  <cpu_time>%f</cpu_time>\n"
	   "</graphics_info>\n",
	   boincv6_skypos_rac,
	   boincv6_skypos_dec,
	   boinc_get_fraction_done(),
	   boinc_worker_thread_cpu_time());
}

int setup_shmem(void) {
  boinc_get_init_data(eah_app_init_data);

  shmem = (HS_SHMEM*)boinc_graphics_make_shmem(EAH_SHMEM_APP_NAME, sizeof(HS_SHMEM));
  if (!shmem) {
    fprintf(stderr, "failed to create shared mem segment\n");
    return(-1);
  }
  update_shmem();
  boinc_register_timer_callback(update_shmem);
  return(0);
}

#endif


#if (BOINC_GRAPHICS == 1) || ((BOINC_GRAPHICS == 2) && defined(_MSC_VER))
int boinc_init_graphics_options(WORKER_FUNC_PTR worker)
{
  return(boinc_init_options_graphics(eah_boinc_options, worker));
}
#endif

void set_boinc_options(void) {
  rcsid = HS_BOINC_OPTIONS_H_RCSID;
  boinc_options_defaults(eah_boinc_options);
#if (BOINC_GRAPHICS > 0)
  // only makes sense on Apps with graphics
  eah_boinc_options.backwards_compatible_graphics = 0;
#endif
}
