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
#include "util.h" /* for BOINC's dtime() */

char* rcsid = "$Id$";
BOINC_OPTIONS eah_boinc_options;
APP_INIT_DATA eah_app_init_data;

/* BOINC APIv6 stuff */

t_progress boincv6_progress;

#ifdef BOINC_APIV6
#include "graphics2.h"

char*shmem = NULL;

static void update_shmem(void) {
  BOINC_STATUS boincstat;

  if (!shmem) return;

  boinc_get_status(&boincstat);

  snprintf(shmem, EAH_SHMEM_SIZE,
	   "<graphics_info>\n"
	   "  <skypos_rac>%f</skypos_rac>\n"
	   "  <skypos_dec>%f</skypos_dec>\n"
	   "  <fraction_done>%f</fraction_done>\n"
	   "  <cpu_time>%f</cpu_time>\n"
	   "  <update_time>%f</update_time>\n"
	   "  <frequency>%f</frequency>\n"
	   "  <bandwidth>%f</bandwidth>\n"
           "  <candidate>\n"
	   "    <frequency>%f<frequency>\n"
	   "    <spindown>%f<spindown>\n"
	   "    <rac>%f<rac>\n"
	   "    <dec>%f<dec>\n"
	   "    <hough_sig>%f<hough_sig>\n"
           "  </candidate>\n"
           "  <boinc_status>\n"
	   "    <no_heartbeat>%d<no_heartbeat>\n"
	   "    <suspended>%d<suspended>\n"
	   "    <quit_request>%d<quit_request>\n"
	   "    <reread_init_data_file>%d<reread_init_data_file>\n"
	   "    <abort_request>%d<abort_request>\n"
	   "    <working_set_size>%f<working_set_size>\n"
	   "    <max_working_set_size>%f<max_working_set_size>\n"
           "  </boinc_status>\n"
	   "</graphics_info>\n",
	   boincv6_progress.skypos_rac,
	   boincv6_progress.skypos_dec,
	   boinc_get_fraction_done(),
	   boinc_worker_thread_cpu_time(),
	   dtime(),
	   boincv6_progress.frequency,
	   boincv6_progress.bandwidth,
	   boincv6_progress.cand_frequency,
	   boincv6_progress.cand_spindown,
	   boincv6_progress.cand_rac,
	   boincv6_progress.cand_dec,
	   boincv6_progress.cand_hough_sign,
	   boincstat.no_heartbeat,
	   boincstat.suspended,
	   boincstat.quit_request,
	   boincstat.reread_init_data_file,
	   boincstat.abort_request,
	   boincstat.working_set_size,
	   boincstat.max_working_set_size);
}

int setup_shmem(void) {
  boinc_get_init_data(eah_app_init_data);

  shmem = (char*)boinc_graphics_make_shmem(EAH_SHMEM_APP_NAME, EAH_SHMEM_SIZE);
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
  boinc_get_init_data(eah_app_init_data);
  boinc_options_defaults(eah_boinc_options);
#if (BOINC_GRAPHICS > 0)
  // only makes sense on Apps with graphics
  eah_boinc_options.backwards_compatible_graphics = 0;
#endif
}
