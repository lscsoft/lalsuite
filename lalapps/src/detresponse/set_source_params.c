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
