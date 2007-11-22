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
#include "graphics_api.h"

char* rcsid = "$Id$";
BOINC_OPTIONS eah_boinc_options;

int boinc_init_graphics_options(WORKER_FUNC_PTR worker)
{
#if (BOINC_GRAPHICS == 1) || ((BOINC_GRAPHICS == 2) && defined(_MSC_VER))
  return(boinc_init_options_graphics(eah_boinc_options, worker));
#endif
}

void set_boinc_options(void) {
  rcsid = HS_BOINC_OPTIONS_H_RCSID;
  boinc_options_defaults(eah_boinc_options);
  eah_boinc_options.backwards_compatible_graphics = 0;
}
