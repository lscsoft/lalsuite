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

#include "hs_boinc_options.h"
#include "boinc_api.h"

char* rcsid = "$Id$";
BOINC_OPTIONS eah_boinc_options;
APP_INIT_DATA eah_app_init_data;

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
  boinc_get_init_data(eah_app_init_data);
  /*
    Just to remind me of the info in there:

    fprintf(f,"User: %s \t\t", init_data.user_name);
    fprintf(f,"Auth: %s \n", init_data.authenticator);
    fprintf(f,"Team: %s \n\n", init_data.team_name);
    fprintf(f,"Project Dir: %s \n", init_data.project_dir);  
    fprintf(f,"BOINC Dir: %s \n", init_data.boinc_dir);  
    fprintf(f,"\n");
    fprintf(f,"Total Credit:\t %f\n", init_data.user_total_credit);
    fprintf(f,"Avg Credit:\t %f\n",     init_data.user_expavg_credit);
    fprintf(f,"\n");
    fprintf(f,"APP Name: %s \n", init_data.app_name);
    fprintf(f,"WU Name:  %s \n", init_data.wu_name);
  */

}
