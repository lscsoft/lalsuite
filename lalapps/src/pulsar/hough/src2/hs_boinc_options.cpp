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

using namespace std;
#include <iostream>

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
  // char buf[1024];

  if (!shmem) return;

  boinc_get_init_data(eah_app_init_data);

  shmem->fraction_done = boinc_get_fraction_done(); // boincv6_fraction_done;
  shmem->skypos_rac    = boincv6_skypos_rac;
  shmem->skypos_dec    = boincv6_skypos_dec;
  shmem->user_credit   = eah_app_init_data.user_total_credit;
  shmem->ravg_credit   = eah_app_init_data.user_expavg_credit;
  shmem->host_credit   = eah_app_init_data.host_total_credit;
  shmem->cpu_time      = boinc_worker_thread_cpu_time();;
  shmem->update_time   = 0; //dtime();
  strncpy(shmem->user_name, eah_app_init_data.user_name, sizeof(shmem->user_name));
  strncpy(shmem->team_name, eah_app_init_data.team_name, sizeof(shmem->team_name));
  strncpy(shmem->app_name,  eah_app_init_data.app_name,  sizeof(shmem->app_name));
  strncpy(shmem->wu_name,   eah_app_init_data.wu_name,   sizeof(shmem->wu_name));
  strncpy(shmem->boincdir,  eah_app_init_data.boinc_dir, sizeof(shmem->boincdir));
  boinc_get_status(&shmem->status);

  /*
    Just to remind me of the info also in init.data

    fprintf(f,"Auth: %s \n", eah_app_init_data.authenticator);
    fprintf(f,"Project Dir: %s \n", eah_app_init_data.project_dir);  
    fprintf(f,"BOINC Dir: %s \n", eah_app_init_data.boinc_dir);  
  */

  char*buf = shmem->xml;

  snprintf(buf, sizeof(buf),
	   "<graphics_info>\n"
	   "  <fraction_done>%f</fraction_done>\n"
	   "  <skypos_rac>%f</skypos_rac>\n"
	   "  <skypos_dec>%f</skypos_dec>\n"
	   "  <user_credit>%f</user_credit>\n"
	   "  <ravg_credit>%f</ravg_credit>\n"
	   "  <host_credit>%f</host_credit>\n"
	   "  <cpu_time>%f</cpu_time>\n"
	   "  <update_time>%f</update_time>\n"
	   "  <user_name>%s</user_name>\n"
	   "  <team_name>%s</team_name>\n"
	   "  <app_name>%s</app_name>\n"
	   "  <wu_name>%s</wu_name>\n"
	   "  <boincdir>%s</boincdir>\n"
	   "</graphics_info>\n",
	   boinc_get_fraction_done(),
	   boincv6_skypos_rac,
	   boincv6_skypos_dec,
	   eah_app_init_data.user_total_credit,
	   eah_app_init_data.user_expavg_credit,
	   eah_app_init_data.host_total_credit,
	   boinc_worker_thread_cpu_time(),
	   0.0, //dtime();
	   eah_app_init_data.user_name,
	   eah_app_init_data.team_name,
	   eah_app_init_data.app_name,
	   eah_app_init_data.wu_name,
	   eah_app_init_data.boinc_dir);
  cout << buf;
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
