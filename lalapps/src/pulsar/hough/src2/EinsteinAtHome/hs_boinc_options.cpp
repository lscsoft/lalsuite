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
#include <boinc/util.h> /* for BOINC's dtime() */
#include <boinc/svn_version.h>

BOINC_OPTIONS eah_boinc_options;
APP_INIT_DATA eah_app_init_data;
int eah_userid, eah_hostid;
char*eah_hostcpid, *eah_username;

/* BOINC APIv6 stuff */

t_progress boincv6_progress;

#ifdef BOINC_APIV6
#include <boinc/graphics2.h>

char*shmem = NULL;

static void update_shmem(void) {
  BOINC_STATUS boincstat;

  if (!shmem) return;

  boinc_get_status(&boincstat);

  memset(shmem, 0, EAH_SHMEM_SIZE);

  int len = snprintf(shmem, EAH_SHMEM_SIZE,
	   "<graphics_info>\n"
	   "  <skypos_rac>%.3f</skypos_rac>\n"
	   "  <skypos_dec>%.3f</skypos_dec>\n"
	   "  <fraction_done>%.3f</fraction_done>\n"
	   "  <cpu_time>%.3f</cpu_time>\n"
	   "  <update_time>%.3f</update_time>\n"
	   "  <frequency>%.3f</frequency>\n"
	   "  <bandwidth>%.3f</bandwidth>\n"
           "  <candidate>\n"
	   "    <frequency>%.3f</frequency>\n"
	   "    <spindown>%.3f</spindown>\n"
	   "    <rac>%.3f</rac>\n"
	   "    <dec>%.3f</dec>\n"
	   "    <hough_sig>%.3f</hough_sig>\n"
           "  </candidate>\n"
           "  <boinc_status>\n"
	   "    <no_heartbeat>%d</no_heartbeat>\n"
	   "    <suspended>%d</suspended>\n"
	   "    <quit_request>%d</quit_request>\n"
	   "    <reread_init_data_file>%d</reread_init_data_file>\n"
	   "    <abort_request>%d</abort_request>\n"
	   "    <working_set_size>%.3f</working_set_size>\n"
	   "    <max_working_set_size>%.3f</max_working_set_size>\n"
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

  if (len >= EAH_SHMEM_SIZE)
    fprintf(stderr, "WARNING: shmem XML hits max size (%d,%u)\n", len, EAH_SHMEM_SIZE);
}

int setup_shmem(void) {
  boinc_get_init_data(eah_app_init_data);

  if (boinc_is_standalone())
    return(0);

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


void set_boinc_options(void) {
  boinc_get_init_data(eah_app_init_data);
#if defined(SVN_REVISION) && (SVN_REVISION >= 22905)
  // this breaks compilation for BOINC older than SVN trunk r22902
  eah_userid   = eah_app_init_data.userid;
#else
  eah_userid   = -1;
#endif
  eah_username = eah_app_init_data.user_name;
  eah_hostid   = eah_app_init_data.hostid;
  eah_hostcpid = eah_app_init_data.host_info.host_cpid;
  boinc_options_defaults(eah_boinc_options);
}
