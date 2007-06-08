/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "set_source_params.h"

void set_source_params(LALSource * source, const char *name, REAL8 ra_rad,
                       REAL8 dec_rad, REAL8 orien_rad)
{
  (void)mystrlcpy(source->name, name, LALNameLength);
  source->equatorialCoords.longitude = ra_rad;
  source->equatorialCoords.latitude  = dec_rad;
  source->equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  source->orientation                = orien_rad;
}
