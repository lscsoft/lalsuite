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

#include "graphics_api.h"

#define HS_BOINC_OPTIONS_H_RCSID "$Id$"

#ifdef __cplusplus
extern "C" {
#endif
  
extern int boinc_init_graphics_options(WORKER_FUNC_PTR);

extern void set_boinc_options(void);

#ifdef __cplusplus
}
#endif
