/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

#ifndef _DETRESPONSE_SETSOURCEPARAMS_H
#define _DETRESPONSE_SETSOURCEPARAMS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/DetResponse.h>

void  set_source_params(LALSource * source, const char * name, REAL8 ra_rad,
                        REAL8 dec_rad, REAL8 orien_rad);

#endif
